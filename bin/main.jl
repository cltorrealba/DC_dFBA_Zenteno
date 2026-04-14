import Pkg

function find_project_root(start::AbstractString)
    d = abspath(start)
    while true
        isfile(joinpath(d, "Project.toml")) && return d
        parent = dirname(d)
        parent == d && return nothing
        d = parent
    end
end

project_root = find_project_root(pwd())
project_root === nothing && error("No se encontro Project.toml desde pwd=$(pwd())")

Pkg.activate(project_root)
Pkg.instantiate()

# Dependencias que pueden faltar en el entorno
for pkg in ["HDF5", "ForwardDiff"]
    if Base.find_package(pkg) === nothing
        println("[env] $(pkg) no encontrado. Instalando...")
        Pkg.add(pkg)
    end
end
Pkg.precompile()

println("[env] Proyecto activo: ", project_root)


# ── Imports ──
using JuMP
using Ipopt
using LinearAlgebra
using DelimitedFiles
using Printf
using Statistics
using Plots
using Libdl
using SparseArrays
using ForwardDiff

# HDF5 para leer semilla
HDF5_AVAILABLE = false
try
    @eval using HDF5
    global HDF5_AVAILABLE = true
    println("[HDF5] cargado OK")
catch e
    @warn "HDF5 no disponible" e
end

# Directorio de salida donde NB1 exporta S, bounds, metadata, semilla
OUT_DIR = get(ENV, "OUT_DIR", joinpath(pwd(), "out"))

S_FILE    = joinpath(OUT_DIR, "S.csv")
LB_FILE   = joinpath(OUT_DIR, "lb.csv")
UB_FILE   = joinpath(OUT_DIR, "ub.csv")
RXN_FILE  = joinpath(OUT_DIR, "rxn_ids.txt")
MET_FILE  = joinpath(OUT_DIR, "met_ids.txt")
META_FILE = joinpath(OUT_DIR, "dfba_vargam_metadata.jl")

for f in [S_FILE, LB_FILE, UB_FILE, RXN_FILE, MET_FILE, META_FILE]
    @assert isfile(f) "Archivo no encontrado: $f"
end
println("[paths] Todos los archivos de NB1 encontrados en: ", OUT_DIR)


# Helpers genéricos (puros y reutilizables)
function _meta_get(meta::Dict{String,Any}, key::String, default)
    haskey(meta, key) ? meta[key] : default
end

function _dict_string_string(x)
    d = Dict{String,String}()
    for (k,v) in pairs(x)
        d[string(k)] = string(v)
    end
    return d
end

function _dict_string_float(x)
    d = Dict{String,Float64}()
    for (k,v) in pairs(x)
        d[string(k)] = Float64(v)
    end
    return d
end

function linear_interp(xq::Vector{Float64}, xp::Vector{Float64}, yp::Vector{Float64})
    isempty(xp) && return zeros(length(xq))
    @assert length(xp) == length(yp)
    yq = similar(xq)
    for i in eachindex(xq)
        x = xq[i]
        if x <= xp[1]
            yq[i] = yp[1]
        elseif x >= xp[end]
            yq[i] = yp[end]
        else
            j = searchsortedlast(xp, x)
            x1, x2 = xp[j], xp[j+1]
            y1, y2 = yp[j], yp[j+1]
            θ = (x - x1) / max(x2 - x1, 1e-12)
            yq[i] = y1 + θ * (y2 - y1)
        end
    end
    return yq
end

function _clip(x, lo, hi)
    return min(max(x, lo), hi)
end

function _finite_diff(values::Vector{Float64}, times::Vector{Float64})
    n = length(values)
    out = zeros(n)
    if n <= 1
        return out
    end
    for k in 1:n
        if k == 1
            dt = max(times[2] - times[1], 1e-9)
            out[k] = (values[2] - values[1]) / dt
        elseif k == n
            dt = max(times[n] - times[n-1], 1e-9)
            out[k] = (values[n] - values[n-1]) / dt
        else
            dt = max(times[k+1] - times[k-1], 1e-9)
            out[k] = (values[k+1] - values[k-1]) / dt
        end
    end
    return out
end

softplus_num(x::Float64, ϵ::Float64) = 0.5 * (x + sqrt(x*x + ϵ*ϵ))
sigmoid_num(x::Float64) = 1.0 / (1.0 + exp(-x))

function _percentiles(x::AbstractVector{<:Real})
    isempty(x) && return Dict("p50"=>0.0, "p90"=>0.0, "p95"=>0.0, "p99"=>0.0)
    xv = sort(Float64.(x))
    q(p) = xv[clamp(ceil(Int, p * length(xv)), 1, length(xv))]
    return Dict("p50"=>q(0.50), "p90"=>q(0.90), "p95"=>q(0.95), "p99"=>q(0.99))
end

function _status(max_res::Real, tol::Real)
    max_res <= tol ? "PASS" : (max_res <= 10tol ? "WARN" : "FAIL")
end

function _family_name(s::Int)
    if s == IDX_X
        return "biomasa"
    elseif s in (IDX_G, IDX_F)
        return "azucares"
    elseif s == IDX_NFREE || s in AA_STATE_IDXS
        return "nitrogeno_y_aa"
    elseif s == IDX_E
        return "etanol"
    elseif s in (IDX_PROT, IDX_CARB)
        return "proteina_carbohidrato"
    elseif s in AROMA_STATE_IDXS
        return "aromas"
    else
        return "otros"
    end
end

function _build_time_maps(hvec::Vector{Float64})
    rad = [0.15505102572168, 0.64494897427832, 1.0]
    t_loc = zeros(length(hvec), 3)
    t0 = 0.0
    for i in eachindex(hvec)
        for j in 1:3
            t_loc[i,j] = t0 + rad[j] * hvec[i]
        end
        t0 += hvec[i]
    end
    return t_loc
end

S = Float64.(readdlm(S_FILE, ','))
lbraw = readdlm(LB_FILE, ',')
ubraw = readdlm(UB_FILE, ',')
RXN_IDS = readlines(RXN_FILE)
MET_IDS = readlines(MET_FILE)

vlb = lbraw isa AbstractVector ? Float64.(lbraw) : Float64.(lbraw[:, 1])
vub = ubraw isa AbstractVector ? Float64.(ubraw) : Float64.(ubraw[:, 1])

nm = size(S, 1)
nv = size(S, 2)

RXN_INDEX = Dict{String, Int}(rid => i for (i, rid) in enumerate(RXN_IDS))
MET_INDEX = Dict{String, Int}(mid => i for (i, mid) in enumerate(MET_IDS))

META = Dict{String,Any}()
if isfile(META_FILE)
    include(META_FILE)
    if @isdefined DFBA_META
        META = deepcopy(DFBA_META)
    else
        # Fallback: old metadata format uses top-level variables, not Dict.
        # Collect critical GAM profile data into META so downstream code works.
        @isdefined(baseline_time_h)       && (META["baseline_time_h"] = baseline_time_h)
        @isdefined(baseline_gam_mmol_gdw) && (META["baseline_gam_mmol_gdw"] = baseline_gam_mmol_gdw)
        println("[META] DFBA_META not found — fallback: collected globals into META")
    end
end

OBJ_ID  = _meta_get(META, "obj_id", "r_2111")
GLU_ID  = _meta_get(META, "glu_id", "r_1714")
FRU_ID  = _meta_get(META, "fru_id", "r_1709")
ETH_ID  = _meta_get(META, "eth_id", "r_1761")
O2_ID   = _meta_get(META, "o2_id",  "r_1992")
ATPM_ID = _meta_get(META, "atpm_id","r_4046")
PROT_RXN_ID = _meta_get(META, "prot_rxn_id", "r_4047")

for rid in [OBJ_ID, GLU_ID, FRU_ID, ETH_ID, O2_ID, ATPM_ID, PROT_RXN_ID]
    @assert haskey(RXN_INDEX, rid) "Falta reacción requerida: $(rid)"
end

obj = RXN_INDEX[OBJ_ID]
glu = RXN_INDEX[GLU_ID]
fru = RXN_INDEX[FRU_ID]
eth = RXN_INDEX[ETH_ID]
o2  = RXN_INDEX[O2_ID]
IDX_ATPM = RXN_INDEX[ATPM_ID]
IDX_PROT_RXN = RXN_INDEX[PROT_RXN_ID]

kinetic_n_default = ["r_1654", "r_1879", "r_1891", "r_1889", "r_1906", "r_1911", "r_1873", "r_1912"]
KINETIC_N_SOURCE_IDS = [string(x) for x in _meta_get(META, "KINETIC_N_SOURCE_IDS", kinetic_n_default)]
for rid in KINETIC_N_SOURCE_IDS
    @assert haskey(RXN_INDEX, rid) "Falta fuente cinética de N: $(rid)"
end
KINETIC_N_SOURCE_IDXS = [RXN_INDEX[rid] for rid in KINETIC_N_SOURCE_IDS]

AA_EXCHANGE_MAP = haskey(META, "aa_exchange_ids") ? _dict_string_string(META["aa_exchange_ids"]) :
    Dict("phe"=>"r_1898", "leu"=>"r_1890", "val"=>"r_1910", "met"=>"r_1893", "tyr"=>"r_1914")
for rid in values(AA_EXCHANGE_MAP)
    @assert haskey(RXN_INDEX, rid) "Falta exchange de AA: $(rid)"
end

AROMA_EXCHANGE_MAP = haskey(META, "aroma_exchange_ids") ? _dict_string_string(META["aroma_exchange_ids"]) :
    Dict("pea"=>"r_1590", "isoamyl"=>"r_1865", "isobutanol"=>"r_1866", "methionol"=>"r_1900", "tyrosol"=>"r_1915")
for rid in values(AROMA_EXCHANGE_MAP)
    @assert haskey(RXN_INDEX, rid) "Falta exchange de aroma/alcohol: $(rid)"
end

AA_KEYS = collect(keys(AA_EXCHANGE_MAP))
AROMA_KEYS = collect(keys(AROMA_EXCHANGE_MAP))
sort!(AA_KEYS)
sort!(AROMA_KEYS)

AA_UPTAKE_IDXS = [RXN_INDEX[AA_EXCHANGE_MAP[k]] for k in AA_KEYS]
AROMA_RXN_IDXS = [RXN_INDEX[AROMA_EXCHANGE_MAP[k]] for k in AROMA_KEYS]

UPTAKE_IDXS = vcat([glu, fru], KINETIC_N_SOURCE_IDXS)
PRODUCT_IDXS = [eth, obj]
n_up = length(UPTAKE_IDXS)
n_prod = length(PRODUCT_IDXS)

IS_GLU = [idx == glu ? 1.0 : 0.0 for idx in UPTAKE_IDXS]
IS_FRU = [idx == fru ? 1.0 : 0.0 for idx in UPTAKE_IDXS]
IS_NIT = [1.0 - IS_GLU[i] - IS_FRU[i] for i in eachindex(UPTAKE_IDXS)]

SELECT_UPTAKE = [Float64(mc == UPTAKE_IDXS[k]) for mc in 1:nv, k in 1:n_up]
SELECT_PRODUCT = [Float64(mc == PRODUCT_IDXS[k]) for mc in 1:nv, k in 1:n_prod]

IS_ETH_prod = [1.0, 0.0]
IS_OBJ_prod = [0.0, 1.0]

N_atoms_map = haskey(META, "n_atoms_map") ? _dict_string_float(META["n_atoms_map"]) : Dict{String,Float64}()
N_frac_map  = haskey(META, "n_frac_map")  ? _dict_string_float(META["n_frac_map"])  : Dict{String,Float64}()

MW_N   = 0.014007
MW_GLU = 0.180156
MW_FRU = 0.180156
MW_ETH = 0.046070
MW_O2  = 0.031998

N_atoms_vec = ones(nv)
N_profile_vec = zeros(nv)
for (rid, val) in pairs(N_atoms_map)
    haskey(RXN_INDEX, rid) && (N_atoms_vec[RXN_INDEX[rid]] = Float64(val))
end
for (rid, val) in pairs(N_frac_map)
    haskey(RXN_INDEX, rid) && (N_profile_vec[RXN_INDEX[rid]] = Float64(val))
end

N_frac = zeros(n_up)
for k in 1:n_up
    idx = UPTAKE_IDXS[k]
    if IS_NIT[k] > 0.5
        N_frac[k] = N_profile_vec[idx] / max(N_atoms_vec[idx] * MW_N, 1e-12)
    end
end

AA_ALPHA_MAP = haskey(META, "aa_alpha") ? _dict_string_float(META["aa_alpha"]) :
    Dict(k => 1.0 / max(length(AA_KEYS), 1) for k in AA_KEYS)
AA_ALPHA_VEC = [get(AA_ALPHA_MAP, k, 0.0) for k in AA_KEYS]

AA_MW_MAP = haskey(META, "aa_mw") ? _dict_string_float(META["aa_mw"]) :
    Dict("phe"=>0.16519, "leu"=>0.13117, "val"=>0.11715, "met"=>0.14921, "tyr"=>0.18119)
AROMA_MW_MAP = haskey(META, "aroma_mw") ? _dict_string_float(META["aroma_mw"]) :
    Dict("pea"=>0.12217, "isoamyl"=>0.08815, "isobutanol"=>0.07412, "methionol"=>0.10619, "tyrosol"=>0.13816)
AROMA_MW_VEC = [get(AROMA_MW_MAP, k, 0.1) for k in AROMA_KEYS]

# Pairwise / ethyl acetate defaults aligned with Notebook 1
pairwise_default = [
    ("r_1862", "r_1865", 0.08),
    ("r_1867", "r_1866", 0.08),
    ("r_2000", "r_1589", 0.08),
]
PAIRWISE_META = haskey(META, "pairwise_constraints") ? META["pairwise_constraints"] : pairwise_default
PAIRWISE_ESTER_IDXS = Int[]
PAIRWISE_ALCOHOL_IDXS = Int[]
PAIRWISE_PHI = Float64[]
for row in PAIRWISE_META
    local ester_rid = string(row[1])
    local alcohol_rid = string(row[2])
    local phi = Float64(row[3])
    if haskey(RXN_INDEX, ester_rid) && haskey(RXN_INDEX, alcohol_rid)
        push!(PAIRWISE_ESTER_IDXS, RXN_INDEX[ester_rid])
        push!(PAIRWISE_ALCOHOL_IDXS, RXN_INDEX[alcohol_rid])
        push!(PAIRWISE_PHI, phi)
    end
end
@assert !isempty(PAIRWISE_ESTER_IDXS) "No se detectaron pares ester↔alcohol."

PAIR_ALCOHOL_SELECT = [Float64(mc == PAIRWISE_ALCOHOL_IDXS[p]) for mc in 1:nv, p in 1:length(PAIRWISE_PHI)]
PAIR_ESTER_SELECT   = [Float64(mc == PAIRWISE_ESTER_IDXS[p])   for mc in 1:nv, p in 1:length(PAIRWISE_PHI)]

EA_SOFT_ESTER_IDX = 0
EA_SOFT_ALCOHOL_IDX = 0
PHI_ETHYL_ACETATE_STATIC = 0.0
if haskey(META, "ethyl_acetate_soft")
    ea = META["ethyl_acetate_soft"]
    ester_rid = string(ea["ester_rid"])
    alcohol_rid = string(ea["alcohol_rid"])
    if haskey(RXN_INDEX, ester_rid) && haskey(RXN_INDEX, alcohol_rid)
        EA_SOFT_ESTER_IDX = RXN_INDEX[ester_rid]
        EA_SOFT_ALCOHOL_IDX = RXN_INDEX[alcohol_rid]
        PHI_ETHYL_ACETATE_STATIC = Float64(ea["phi"])
    end
else
    if haskey(RXN_INDEX, "r_1765") && haskey(RXN_INDEX, ETH_ID)
        EA_SOFT_ESTER_IDX = RXN_INDEX["r_1765"]
        EA_SOFT_ALCOHOL_IDX = RXN_INDEX[ETH_ID]
        PHI_ETHYL_ACETATE_STATIC = 0.005
    end
end

EA_ALCOHOL_SELECT = [Float64(mc == EA_SOFT_ALCOHOL_IDX) for mc in 1:nv]
EA_ESTER_SELECT   = [Float64(mc == EA_SOFT_ESTER_IDX) for mc in 1:nv]
OBJ_SELECT  = [Float64(mc == obj) for mc in 1:nv]
ATPM_SELECT = [Float64(mc == IDX_ATPM) for mc in 1:nv]

println("nm=$(nm), nv=$(nv), nAA=$(length(AA_KEYS)), nAroma=$(length(AROMA_KEYS)), nPair=$(length(PAIRWISE_PHI))")
println("[META] keys loaded: $(length(META)) entries")

# -----------------------------
# Parámetros cinéticos / composición
# Homologados a NB1 (paper2010)
# -----------------------------
MU0_nom    = 0.18
YXN_nom    = 19.69
YXG_nom    = 1.60
YXF_nom    = 1.60
YEG_nom    = 0.49
YEF_nom    = 0.49
Kn0_nom    = 0.01
Kg0_nom    = 7.5
Kf0_nom    = 7.5
Kig0_nom   = 55.0
Kie0_nom   = 40.0
Kd0_nom    = 0.00044
betaG0_nom = 0.225
betaF0_nom = 0.225
MRATE_0    = 0.01

MU0 = MU0_nom
YEG = YEG_nom
YEF = YEF_nom
YXN = YXN_nom
R = 8.314
EPS = 1e-9

PROT_CONTENT_0 = Float64(_meta_get(META, "PROT_CONTENT_0", 0.46))
CARB_CONTENT_0 = Float64(_meta_get(META, "CARB_CONTENT_0", 0.37))
RNA_FRAC       = Float64(_meta_get(META, "RNA_FRAC", 0.06))
K_DEATH        = Float64(_meta_get(META, "K_DEATH", 0.005))
TURNOVER_LAMBDA = Float64(_meta_get(META, "TURNOVER_LAMBDA", 0.03))
XA_FRACTION     = Float64(_meta_get(META, "XA_FRACTION", 1.0))
N_AMMONIA_FRACTION = Float64(_meta_get(META, "N_AMMONIA_FRACTION", 0.50))
N_TOTAL_DEPLETION_THRESHOLD = Float64(_meta_get(META, "N_TOTAL_DEPLETION_THRESHOLD", 1e-3))
K_AA_UPTAKE_GROWTH = Float64(_meta_get(META, "K_AA_UPTAKE_GROWTH", 0.08))

ATPM_LB_NO_GROWTH = Float64(_meta_get(META, "ATPM_LB_NO_GROWTH", 0.70))
ATPM_UB_NO_GROWTH = Float64(_meta_get(META, "ATPM_UB_NO_GROWTH", 1000.0))

GAM_BASE = Float64(_meta_get(META, "GAM_BASE", 24.7))
GAM_COEFF_P = Float64(_meta_get(META, "GAM_COEFF_P", 16.965))
GAM_COEFF_R = Float64(_meta_get(META, "GAM_COEFF_R", 1.638))
GAM_COEFF_C = Float64(_meta_get(META, "GAM_COEFF_C", 5.210))
Pbase_global = Float64(_meta_get(META, "Pbase_global", PROT_CONTENT_0))
Cbase_global = Float64(_meta_get(META, "Cbase_global", CARB_CONTENT_0))
Rbase_global = Float64(_meta_get(META, "Rbase_global", RNA_FRAC))

function compute_full_gam(P, Rna, Carb, Pbase, Rbase, Cbase)
    Pfactor = P / max(Pbase, 1e-9)
    Rfactor = Rna / max(Rbase, 1e-9)
    Cfactor = max(0.0, (Cbase + Pbase - P - Rna) / max(Cbase, 1e-9))
    return GAM_BASE + GAM_COEFF_P * Pfactor + GAM_COEFF_R * Rfactor + GAM_COEFF_C * Cfactor
end

GAM_REF = compute_full_gam(PROT_CONTENT_0, RNA_FRAC, CARB_CONTENT_0, Pbase_global, Rbase_global, Cbase_global)

# -----------------------------
# Perfil térmico e inyección
# -----------------------------
T_BASE  = try parse(Float64, get(ENV, "T_CONST", "293.15")) catch; 293.15 end
T_STEPS = [36.0, 96.0]
T_DELTAS = [5.0, 3.0]
T_STEEP = 0.5

function dynamic_temperature(t)
    val = T_BASE
    for i in eachindex(T_STEPS)
        σ = 1.0 / (1.0 + exp(-T_STEEP * (t - T_STEPS[i])))
        val += T_DELTAS[i] * σ
    end
    return val
end

function death_rate_T(E, T_val)
    Td = -0.0001 * E^3 + 0.0049 * E^2 - 0.1279 * E + 315.89
    s = 0.5 * (1.0 + tanh(0.5 * (T_val - Td)))
    base = Kd0_nom * exp(0.0415 * E + (130000.0 * (T_val - 305.65)) / (305.65 * R * T_val))
    return base * s
end

SQRT_2PI = sqrt(2.0 * pi)
function smooth_injection(t, t_shot, dose, width)
    # Sin early-return con comparacion: incompatible con ForwardDiff Dual numbers.
    # La gaussiana decae a ~1e-11 a 5*width, numéricamente equivale a 0.
    return (dose / (width * SQRT_2PI)) * exp(-0.5 * ((t - t_shot) / width)^2)
end

T_INJ_1 = 0.0; DOSE_1 = 0.0; WIDTH_1 = 5.0
T_INJ_2 = 0.0; DOSE_2 = 0.0; WIDTH_2 = 5.0

# -----------------------------
# GAM exógena desde Notebook 1 (defer hasta tener malla nfe/th)
# -----------------------------

function _vector_from_meta(meta::Dict{String,Any}, key::String)
    if !haskey(meta, key)
        return Float64[]
    end
    return [Float64(x) for x in meta[key]]
end

baseline_time_h = _vector_from_meta(META, "baseline_time_h")
baseline_gam_mmol_gdw = _vector_from_meta(META, "baseline_gam_mmol_gdw")

function build_gam_profiles(nfe_local::Int, th_local::Float64, gam_ref::Float64)
    h_local = fill(th_local / nfe_local, nfe_local)
    tfe = cumsum(h_local)
    if !isempty(baseline_time_h) && !isempty(baseline_gam_mmol_gdw)
        gam_fe = linear_interp(Float64.(tfe), baseline_time_h, baseline_gam_mmol_gdw)
    else
        gam_fe = fill(gam_ref, nfe_local)
    end
    gam_extra = max.(0.0, gam_fe .- gam_ref)
    return gam_fe, gam_extra
end

println(@sprintf("[GAM] baseline points: time=%d gam=%d", length(baseline_time_h), length(baseline_gam_mmol_gdw)))
println("[GAM] Perfil se evaluará en sección de construcción de modelo (nfe/th canónicos).")

# Configuracion COIN-HSL identica a sparsepatch
const HSL_BIN_DIR = raw"C:\COIN_HSL\CoinHSL.v2023.11.17.x86_64-w64-mingw32-libgfortran5\bin"
const HSL_DLL_FILE = joinpath(HSL_BIN_DIR, "libcoinhsl.dll")

ENV["IPOPT_HSLLIB"] = HSL_DLL_FILE
ENV["HSL_DLL_PATH"] = HSL_DLL_FILE

path_sep = Sys.iswindows() ? ';' : ':'
path_entries = split(get(ENV, "PATH", ""), path_sep; keepempty=false)
if !(HSL_BIN_DIR in path_entries)
    ENV["PATH"] = HSL_BIN_DIR * path_sep * get(ENV, "PATH", "")
end

if isfile(HSL_DLL_FILE)
    try
        h = Libdl.dlopen(HSL_DLL_FILE); Libdl.dlclose(h)
        println("[HSL] libcoinhsl.dll OK")
    catch err
        @warn "HSL dll no se pudo abrir" err
    end
else
    @warn "HSL_DLL_FILE no existe" HSL_DLL_FILE
end
println("[HSL] hsllib = ", ENV["IPOPT_HSLLIB"])

# ─────────────────────────────────────────────────────────────────────────
# Malla temporal (Radau IIA, 3 puntos de colocacion por elemento finito)
# ─────────────────────────────────────────────────────────────────────────
nfe = try parse(Int, get(ENV, "NFE", "18")) catch; 18 end
ncp = 3
th  = try parse(Float64, get(ENV, "TH", "72.0")) catch; 168.0 end
h   = th / nfe
hm  = fill(h, nfe)
var_h = 0.50

# ─────────────────────────────────────────────────────────────────────────
# Indices de estado (deben coincidir con el .jl v2 y con celda 46 de NB1)
# ─────────────────────────────────────────────────────────────────────────
IDX_X     = 1
IDX_NFREE = 2
IDX_G     = 3
IDX_F     = 4
IDX_E     = 5
IDX_O2    = 6
IDX_PROT  = 7
IDX_CARB  = 8

AA_STATE_IDXS = collect(9:(8 + length(AA_KEYS)))
AROMA_STATE_IDXS = collect((9 + length(AA_KEYS)):(8 + length(AA_KEYS) + length(AROMA_KEYS)))

# Estado adicional v2: pool de N reciclado desde turnover de proteina
IDX_NREC = 8 + length(AA_KEYS) + length(AROMA_KEYS) + 1
nc = IDX_NREC   # total de estados

println("nc = $nc  (8 base + $(length(AA_KEYS)) AA + $(length(AROMA_KEYS)) aroma + 1 N_rec)")

# ─────────────────────────────────────────────────────────────────────────
# Constantes v2 (pool de N reciclado y ATPM en fase crecimiento)
# ─────────────────────────────────────────────────────────────────────────
Y_N_FROM_PROT      = 0.16    # fraccion de N liberado por turnover de proteina
K_NREC_UPTAKE      = 0.05    # tasa de drenaje de N_rec para fuelear uptake AA en turnover
BAL_ATPM_GROWTH_LB = 0.7     # piso ATPM durante fase de crecimiento
BAL_ATPM_GROWTH_UB = 1000.0  # techo ATPM durante fase de crecimiento

# ─────────────────────────────────────────────────────────────────────────
# Condiciones iniciales
# ─────────────────────────────────────────────────────────────────────────
X0 = 0.5
N0_total = 0.14
N0_ammonia = N0_total * N_AMMONIA_FRACTION
N0_from_aa = max(0.0, N0_total - N0_ammonia)
AA0_each = (N0_from_aa / MW_N) / max(length(AA_KEYS), 1)

c0 = zeros(nc)
c0[IDX_X] = X0
c0[IDX_NFREE] = N0_ammonia
c0[IDX_G] = 110.0
c0[IDX_F] = 110.0
c0[IDX_E] = 0.0
c0[IDX_O2] = 0.0
c0[IDX_PROT] = X0 * PROT_CONTENT_0
c0[IDX_CARB] = X0 * CARB_CONTENT_0
for a in eachindex(AA_STATE_IDXS)
    c0[AA_STATE_IDXS[a]] = AA0_each
end
for a in eachindex(AROMA_STATE_IDXS)
    c0[AROMA_STATE_IDXS[a]] = 0.0
end
c0[IDX_NREC] = 0.0  # pool N reciclado inicia vacio

println("c0: ", c0)

# ─────────────────────────────────────────────────────────────────────────
# Escalas de estado (cs) -- normalizacion para mejorar condicionamiento
# ─────────────────────────────────────────────────────────────────────────
cs = ones(nc)
cs[IDX_NFREE] = 0.2
cs[IDX_G] = 100.0
cs[IDX_F] = 100.0
cs[IDX_E] = 10.0
cs[IDX_O2] = 0.01
cs[IDX_PROT] = max(0.1, c0[IDX_PROT])
cs[IDX_CARB] = max(0.1, c0[IDX_CARB])
for s in AA_STATE_IDXS
    cs[s] = max(0.1, c0[s])
end
for s in AROMA_STATE_IDXS
    cs[s] = 0.01
end
cs[IDX_NREC] = 0.01  # escala N reciclado (debe coincidir con NB1 cell 48)

# ─────────────────────────────────────────────────────────────────────────
# Escalas de flujos (vs) -- reescala flujos con bounds grandes
# ─────────────────────────────────────────────────────────────────────────
FLUX_SCALE_TARGET = 50.0
vs = ones(nv)
for rx in 1:nv
    br = max(abs(vlb[rx]), abs(vub[rx]))
    if br > FLUX_SCALE_TARGET
        vs[rx] = br / FLUX_SCALE_TARGET
    end
end

# ─────────────────────────────────────────────────────────────────────────
# Selectores KKT (matrices/vectores auxiliares para stationarity)
# ─────────────────────────────────────────────────────────────────────────
SELECT_AAUPTAKE = [Float64(mc == AA_UPTAKE_IDXS[a]) for mc in 1:nv, a in eachindex(AA_UPTAKE_IDXS)]
D_GROWTH = zeros(nv); D_GROWTH[obj] = -1.0
D_TURNOVER = zeros(nv); D_TURNOVER[IDX_ATPM] = -1.0
D_AAUP = zeros(nv)
for idx in AA_UPTAKE_IDXS
    D_AAUP[idx] = 1e-3
end

# ─────────────────────────────────────────────────────────────────────────
# Perfil GAM exogeno sobre la malla de elementos finitos
# ─────────────────────────────────────────────────────────────────────────
GAM_FE, GAM_EXTRA_FE = build_gam_profiles(nfe, th, GAM_REF)

# ─────────────────────────────────────────────────────────────────────────
# Registro opcional del paquete HSL_jll licenciado
# ─────────────────────────────────────────────────────────────────────────
HSL_JLL_PATH = strip(get(ENV, "HSL_JLL_PATH", ""))
if !isempty(HSL_JLL_PATH)
    try
        Pkg.develop(path = HSL_JLL_PATH)
        println("HSL_jll registrado desde HSL_JLL_PATH = ", HSL_JLL_PATH)
    catch err
        @warn "No se pudo registrar HSL_jll desde HSL_JLL_PATH." exception = (err, catch_backtrace())
    end
end

# ─────────────────────────────────────────────────────────────────────────
# Parametros IPOPT / MPCC
# ─────────────────────────────────────────────────────────────────────────
IPOPT_LINEAR_SOLVER = "ma86"
IPOPT_PRINT_LEVEL = try parse(Int, get(ENV, "IPOPT_PRINT_LEVEL", "5")) catch; 5 end
IPOPT_TOL = try parse(Float64, get(ENV, "IPOPT_TOL", "1e-4")) catch; 1e-4 end
IPOPT_ACCEPTABLE_TOL = try parse(Float64, get(ENV, "IPOPT_ACCEPTABLE_TOL", "1e-2")) catch; 1e-2 end
IPOPT_ACCEPTABLE_ITER = try parse(Int, get(ENV, "IPOPT_ACCEPTABLE_ITER", "12")) catch; 12 end
# Homotopia Scholtes: cada nivel de eps re-resuelve con warm start. Max iters por nivel.
IPOPT_MAX_ITER = try parse(Int, get(ENV, "IPOPT_MAX_ITER", "300")) catch; 300 end
IPOPT_CONSTR_VIOL_TOL = try parse(Float64, get(ENV, "IPOPT_CONSTR_VIOL_TOL", "1e-5")) catch; 1e-5 end
IPOPT_COMPL_INF_TOL = try parse(Float64, get(ENV, "IPOPT_COMPL_INF_TOL", "1e-4")) catch; 1e-4 end
IPOPT_MUMPS_MEM_PERCENT = try parse(Int, get(ENV, "IPOPT_MUMPS_MEM_PERCENT", "20")) catch; 20 end

# IPOPT no soporta rutas con caracteres no-ASCII (ñ, ó, etc.)
# Usamos un directorio temporal con ruta ASCII para el log.
_ipopt_tmp = joinpath(tempdir(), "ipopt_dfba")
mkpath(_ipopt_tmp)
IPOPT_OUTPUT_FILE = joinpath(_ipopt_tmp, "ipopt_log_main.txt")
println("[IPOPT] output_file = ", IPOPT_OUTPUT_FILE)

# Pesos de penalizacion complementariedad (LEGACY, NO USADOS con Scholtes).
# Las complementariedades ahora se manejan via bounds |FO_*| <= eps y homotopia.
PHI_L = 0.0
PHI_U = 0.0
PHI_UPT = 0.0
PHI_PROD = 0.0
PHI_AA = 0.0
PHI_PAIR = 0.0
PHI_EA = 0.0
PHI_ATPM = 0.0
# Regularizacion outer pFBA: Q_REG * ||v_phys||^2. Debe ser pequena pero no 0,
# para dar gradiente informativo a IPOPT.
Q_REG = 1e-4
FLUX_SMOOTH_WEIGHT = 0.0
SOFTPLUS_V_EPS = 1e-6
PHASE_SMOOTH_EPS = 5e-4

# Schedule homotopia Scholtes: |FO_*| <= eps_k decreciente
EPS_SCHEDULE = [1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5]

APPLY_PRODUCT_CAPS = lowercase(get(ENV, "APPLY_PRODUCT_CAPS", "false")) in ("1", "true", "yes", "on")
EPS_FLUX = try parse(Float64, get(ENV, "EPS_FLUX", "1e-5")) catch; 1e-5 end

println("IPOPT linear solver = ", IPOPT_LINEAR_SOLVER)
println("nc=$(nc), nfe=$(nfe), ncp=$(ncp), th=$(th), apply_caps=$(APPLY_PRODUCT_CAPS)")
println(@sprintf("GAM_FE range = [%.4f, %.4f] | GAM_EXTRA range = [%.4f, %.4f]",
        minimum(GAM_FE), maximum(GAM_FE), minimum(GAM_EXTRA_FE), maximum(GAM_EXTRA_FE)))

# Lectura de seed_mpcc_from_cell46.h5 producido por NB1.
SEED_FILE = joinpath(OUT_DIR, "seed_mpcc_from_cell46.h5")

function _load_seed_h5(path::String)
    isfile(path) || (@warn "SEED_FILE no existe" path; return nothing)
    if !HDF5_AVAILABLE
        @warn "HDF5.jl no cargado; instala/activa entorno y vuelve a correr celda 1 y 2."
        return nothing
    end
    d = Dict{Symbol,Any}()
    h5open(path, "r") do f
        for k in keys(f)
            d[Symbol(k)] = read(f[k])
        end
    end
    return d
end

function _to_tuple(x)
    x isa Tuple && return x
    return Tuple(x)
end

function _coerce_seed_shape(name::Symbol, arr, expected)
    expected_t = _to_tuple(expected)
    sz = size(arr)

    # Exact match
    if sz == expected_t
        return arr, false
    end

    # 2D transpuesta (h5py <-> Julia)
    if length(expected_t) == 2 && ndims(arr) == 2
        if sz == (expected_t[2], expected_t[1])
            return permutedims(arr, (2, 1)), true
        end
    end

    # 3D con ejes invertidos (ej: (ncp,nfe,nc) -> (nc,nfe,ncp))
    if length(expected_t) == 3 && ndims(arr) == 3
        if sz == (expected_t[3], expected_t[2], expected_t[1])
            return permutedims(arr, (3, 2, 1)), true
        end
    end

    # Vector almacenado como matriz fila/columna
    if length(expected_t) == 1 && ndims(arr) == 2
        if sz == (1, expected_t[1]) || sz == (expected_t[1], 1)
            return vec(arr), true
        end
    end

    return arr, false
end

seed_raw = _load_seed_h5(SEED_FILE)

# Dimensiones esperadas por el modelo Julia (para corregir ejes si vienen de h5py)
expected_shapes = Dict{Symbol,Any}(
    :c => (nc, nfe, ncp),
    :cdot => (nc, nfe, ncp),
    :v => (nv, nfe),
    :hv => (nfe,),
    :lambda_ => (nm, nfe),
    :alpha_L => (nv, nfe),
    :alpha_U => (nv, nfe),
    :alpha_upt => (n_up, nfe),
    :alpha_prod => (n_prod, nfe),
    :alpha_aa => (length(AA_KEYS), nfe),
    :alpha_pair => (length(PAIRWISE_PHI), nfe),
    :alpha_ea => (nfe,),
    :alpha_atpm_floor => (nfe,),
    :FO_L => (nv, nfe),
    :FO_U => (nv, nfe),
    :FO_upt => (n_up, nfe),
    :FO_prod => (n_prod, nfe),
    :FO_aa => (length(AA_KEYS), nfe),
    :FO_pair => (length(PAIRWISE_PHI), nfe),
    :FO_ea => (nfe,),
    :FO_atpm => (nfe,),
)

# ── Normalizar orientacion + reescalar semilla (si corresponde) ──
if seed_raw !== nothing
    seed = Dict{Symbol,Any}()

    for (k, v) in seed_raw
        vv = copy(v)
        if haskey(expected_shapes, k)
            vv2, changed = _coerce_seed_shape(k, vv, expected_shapes[k])
            if changed
                println("[seed] ejes corregidos en $(k): $(size(vv)) -> $(size(vv2))")
            end
            vv = vv2
        end
        seed[k] = vv
    end

    # ── NB1 cell 48 ya exporta c, cdot, v en unidades ESCALADAS (c/cs, v/vs).
    #    NO re-escalar aqui; hacerlo causaba doble-escalamiento y inf_pr = 1e4.
    println("[seed] semilla ya en unidades escaladas (NB1 cell 48) → sin re-escalar")

    println("[seed] claves cargadas (post-normalizacion):")
    for (k, v) in seed
        println("  ", rpad(String(k), 20), " shape=", size(v))
    end
else
    seed = nothing
    println("[seed] no disponible - se corre sin warm-start de semilla.")
end

# Preview rapido: head de cada variable de la semilla
if seed !== nothing
    for k in (:c, :cdot, :v, :hv, :lambda_, :alpha_L, :alpha_U,
              :alpha_upt, :alpha_prod, :alpha_aa, :alpha_pair,
              :FO_L, :FO_U, :FO_upt, :FO_prod, :FO_aa, :FO_pair,
              :FO_ea, :FO_atpm, :alpha_ea, :alpha_atpm_floor)
        if haskey(seed, k)
            A = seed[k]
            if length(A) > 0
                println("-- $k  shape=$(size(A))  range=[$(minimum(A)), $(maximum(A))]")
            else
                println("-- $k  shape=$(size(A))  (vacio)")
            end
        end
    end
    # Mostrar claves en seed que NO estan en la lista esperada
    expected = Set([:c, :cdot, :v, :hv, :lambda_, :alpha_L, :alpha_U,
                    :alpha_upt, :alpha_prod, :alpha_aa, :alpha_pair,
                    :FO_L, :FO_U, :FO_upt, :FO_prod, :FO_aa, :FO_pair,
                    :FO_ea, :FO_atpm, :alpha_ea, :alpha_atpm_floor])
    extras = setdiff(Set(keys(seed)), expected)
    if !isempty(extras)
        println("
Claves extra en seed (no esperadas por MPCC): ", extras)
    end
else
    println("Sin semilla cargada.")
end


# ─────────────────────────────────────────────────────────────────────────
# Fix de factibilidad primal del seed:
#   (A) Proyecta v al nulo de S (S·v_phys = 0) via min-norm correction
#   (B) Clipa v a [vlb, vub] (en unidades fisicas)
#   (C) Recomputa cdot desde c satisfaciendo colocacion Radau exactamente
# Debe llamarse ANTES de reconstruct_duals_from_seed!
# ─────────────────────────────────────────────────────────────────────────
function fix_seed_feasibility!(seed::Dict)
    if !(haskey(seed, :v) && haskey(seed, :c) && haskey(seed, :cdot) && haskey(seed, :hv))
        @warn "[fix-seed] faltan claves; skip"
        return seed
    end
    v_s    = seed[:v]
    c_s    = seed[:c]
    cdot_s = seed[:cdot]
    hv_s   = seed[:hv]

    # ---- (A) Proyeccion de v al nulo de S ----
    # S de yeast-GEM es rank-deficiente (nm=2806, rank < nm), asi que (S·diag(vs²)·S')
    # es singular. Se usa eigendecomposicion truncada para obtener la pseudoinversa
    # genuina en vez de Tikhonov (que undershoots la correccion).
    S_mat = sparse(S)
    vs_diag = Diagonal(vs)
    S_eff_T = S_mat * vs_diag             # (nm, nv)
    SeSe_dense = Matrix(S_eff_T * S_eff_T')  # (nm, nm) — para eigen
    println("[fix-seed] (A) Computando eigendecomposicion de S·diag(vs²)·S' ($(nm)×$(nm))...")
    E = eigen(Symmetric(SeSe_dense))
    eig_thresh = maximum(abs.(E.values)) * 1e-10
    n_rank = count(v -> abs(v) > eig_thresh, E.values)
    pinv_vals = [abs(v) > eig_thresh ? 1.0/v : 0.0 for v in E.values]
    SeSe_pinv = E.vectors * Diagonal(pinv_vals) * E.vectors'
    println("[fix-seed] (A) rank(S·diag(vs²)·S') = $(n_rank) / $(nm)  ($(nm - n_rank) modos nulos truncados)")

    # ---- (A+B) QP per FE: min ||v - v_seed||²  s.t. S_indep·v_phys=0, vlb≤v_phys≤vub ----
    # S has rank 2593/2806 — using all 2806 rows makes the QP KKT system singular,
    # causing IPOPT to hit ITERATION_LIMIT. Solution: project S onto its independent
    # rows using the eigendecomposition we already computed.
    # V_r' * S has full row rank = n_rank, and null(V_r'*S) = null(S).
    eig_mask = [abs(v) > eig_thresh for v in E.values]
    V_r = E.vectors[:, eig_mask]          # (nm × n_rank)
    S_mat_dense = Matrix(S_mat)            # (nm × nv)
    S_indep = V_r' * S_mat_dense           # (n_rank × nv), full row rank
    println("[fix-seed] (A+B) S_indep: $(size(S_indep)) ($(n_rank) independent constraints)")

    sc_before_qp = 0.0
    sc_after_qp  = 0.0
    for i in 1:nfe
        v_seed_phys = [v_s[rx,i] * vs[rx] for rx in 1:nv]
        sc_before_qp = max(sc_before_qp, norm(S_mat * v_seed_phys, Inf))

        qp = Model(Ipopt.Optimizer)
        set_silent(qp)
        set_optimizer_attribute(qp, "linear_solver", IPOPT_LINEAR_SOLVER)
        set_optimizer_attribute(qp, "hsllib", get(ENV, "IPOPT_HSLLIB", ""))
        set_optimizer_attribute(qp, "max_iter", 500)
        set_optimizer_attribute(qp, "tol", 1e-8)
        set_optimizer_attribute(qp, "constr_viol_tol", 1e-10)
        @variable(qp, vlb[rx] <= vp[rx=1:nv] <= vub[rx], start = clamp(v_seed_phys[rx], vlb[rx], vub[rx]))
        @constraint(qp, S_indep * vp .== 0)
        @objective(qp, Min, sum((vp[rx] - v_seed_phys[rx])^2 for rx in 1:nv))
        optimize!(qp)

        st = termination_status(qp)
        if st in (MOI.LOCALLY_SOLVED, MOI.OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED,
                  MOI.ALMOST_OPTIMAL)
            vp_sol = value.(vp)
            for rx in 1:nv
                v_s[rx,i] = vp_sol[rx] / vs[rx]
            end
            sc_after_qp = max(sc_after_qp, norm(S_mat * vp_sol, Inf))
        else
            @warn "[fix-seed] QP FE=$i status=$st; keeping original v"
            sc_after_qp = max(sc_after_qp, norm(S_mat * v_seed_phys, Inf))
        end
        if i == 1 || i == nfe || i % 5 == 0
            println("[fix-seed] (A+B) QP FE $i/$(nfe) done: status=$st")
        end
    end
    println("[fix-seed] (A+B) QP projection complete." *
            " max|Sv| before=$(round(sc_before_qp, sigdigits=4))" *
            " after=$(round(sc_after_qp, sigdigits=4))")

    # ---- (C) Recomputar cdot para satisfacer colocacion Radau exactamente ----
    radau_colmat = [0.19681547722366  -0.06553542585020   0.02377097434822;
                    0.39442431473909   0.29207341166523  -0.04154875212600;
                    0.37640306270047   0.51248582618842   0.11111111111111]
    colmat_inv = inv(radau_colmat)

    c0s_local = [c0[l] / cs[l] for l in 1:nc]
    coll_max_before = 0.0
    coll_max_after  = 0.0
    for l in 1:nc, i in 1:nfe
        c_prev = (i == 1) ? c0s_local[l] : c_s[l, i-1, ncp]
        delta_c = [c_s[l,i,j] - c_prev for j in 1:ncp]
        for j in 1:ncp
            rhs = c_prev
            for k in 1:ncp
                rhs += hv_s[i] * radau_colmat[j,k] * cdot_s[l,i,k]
            end
            coll_max_before = max(coll_max_before, abs(c_s[l,i,j] - rhs))
        end
        cdot_new = colmat_inv * delta_c ./ hv_s[i]
        for k in 1:ncp
            cdot_s[l,i,k] = cdot_new[k]
        end
    end
    for l in 1:nc, i in 1:nfe
        c_prev = (i == 1) ? c0s_local[l] : c_s[l, i-1, ncp]
        for j in 1:ncp
            rhs = c_prev
            for k in 1:ncp
                rhs += hv_s[i] * radau_colmat[j,k] * cdot_s[l,i,k]
            end
            coll_max_after = max(coll_max_after, abs(c_s[l,i,j] - rhs))
        end
    end
    println("[fix-seed] (C) cdot recomputed: coll resid before=$(round(coll_max_before, sigdigits=4))" *
            " after=$(round(coll_max_after, sigdigits=4))")

    seed[:v]    = v_s
    seed[:cdot] = cdot_s
    return seed
end

if seed !== nothing
    fix_seed_feasibility!(seed)
end

# ─────────────────────────────────────────────────────────────────────────
# Reconstruccion de duales (opcion b): dado (c, v) de la semilla, resuelve
# por mínimos cuadrados la estacionariedad Lagrangiana por cada FE asumiendo
# alpha_* = 0 y ajustando lambda_. Luego una pasada de active-set asigna
# alpha_L/alpha_U donde v esta cerca de sus cotas, reduciendo el residuo
# inicial de stationarity antes de entrar a IPOPT.
# ─────────────────────────────────────────────────────────────────────────
function reconstruct_duals_from_seed!(seed::Dict)
    if !(haskey(seed, :v) && haskey(seed, :c))
        @warn "[dual-seed] seed sin :v o :c; skip reconstruccion"
        return seed
    end
    v_scaled = seed[:v]
    c_scaled = seed[:c]
    if size(v_scaled) != (nv, nfe) || size(c_scaled) != (nc, nfe, ncp)
        @warn "[dual-seed] shapes inesperados" v=size(v_scaled) c=size(c_scaled)
        return seed
    end

    # ---- phase_g_fe a partir de la semilla (punto ncp de cada FE) ----
    phase_g_fe_v = zeros(nfe)
    for i in 1:nfe
        cNfree_val = c_scaled[IDX_NFREE, i, ncp] * cs[IDX_NFREE]
        cNaa_val = 0.0
        for a in 1:length(AA_STATE_IDXS)
            cNaa_val += c_scaled[AA_STATE_IDXS[a], i, ncp] * cs[AA_STATE_IDXS[a]]
        end
        cNaa_val *= MW_N
        cNtot_val = cNfree_val + cNaa_val
        phase_g_fe_v[i] = 1.0 / (1.0 + exp(-(cNtot_val - N_TOTAL_DEPLETION_THRESHOLD) / PHASE_SMOOTH_EPS))
    end
    phase_t_fe_v = 1.0 .- phase_g_fe_v

    # ---- d_phase[mc,i] ----
    d_phase_v = zeros(nv, nfe)
    for mc in 1:nv, i in 1:nfe
        d_phase_v[mc,i] = phase_g_fe_v[i] * D_GROWTH[mc] +
                          phase_t_fe_v[i] * D_TURNOVER[mc] +
                          phase_g_fe_v[i] * D_AAUP[mc]
    end

    # ---- LSQ: S^T lambda = -(d_phase + Q_REG*v_phys) por FE ----
    S_mat = Matrix(S)          # (nm, nv)
    ST = Matrix(S_mat')        # (nv, nm)
    lambda_new = zeros(nm, nfe)
    alpha_L_new = zeros(nv, nfe)
    alpha_U_new = zeros(nv, nfe)
    res_before = Float64[]
    res_after  = Float64[]
    tol_bound = 1e-4
    for i in 1:nfe
        b = zeros(nv)
        for mc in 1:nv
            v_phys_mc = v_scaled[mc, i] * vs[mc]
            b[mc] = -(d_phase_v[mc,i] + Q_REG * v_phys_mc)
        end
        push!(res_before, norm(b))
        lam_i = ST \ b
        lambda_new[:, i] .= lam_i

        resid = ST * lam_i - b  # stationarity residual con alpha=0
        # active-set: asignar alpha_L/alpha_U donde v esta cerca de cota y resid != 0
        for mc in 1:nv
            v_phys_mc = v_scaled[mc, i] * vs[mc]
            near_lb = (v_phys_mc - vlb[mc]) < tol_bound * max(abs(vlb[mc]), 1.0)
            near_ub = (vub[mc] - v_phys_mc) < tol_bound * max(abs(vub[mc]), 1.0)
            r = resid[mc]
            # stationarity: ... + alpha_L + alpha_U + ... = 0
            #   => alpha_L + alpha_U = -r
            if near_ub && !near_lb && r < 0
                alpha_U_new[mc,i] = -r          # alpha_U >= 0
            elseif near_lb && !near_ub && r > 0
                alpha_L_new[mc,i] = -r          # alpha_L <= 0
            end
        end

        # residuo final post-activeset
        r_final = copy(resid)
        for mc in 1:nv
            r_final[mc] += alpha_L_new[mc,i] + alpha_U_new[mc,i]
        end
        push!(res_after, norm(r_final))
    end

    seed[:lambda_] = lambda_new
    seed[:alpha_L] = alpha_L_new
    seed[:alpha_U] = alpha_U_new

    # ---- recomputar FO_L/FO_U/FO_* consistentes con los NUEVOS multiplicadores ----
    # CRITICO: si dejamos los FO_* viejos del seed, las restricciones de igualdad
    #   FO_U == (v_phys - vub) * alpha_U
    # quedan violadas masivamente (inf_pr explota). Debemos recomputarlos desde
    # los (v, alpha) que vamos a sembrar.
    FO_L_new = zeros(nv, nfe)
    FO_U_new = zeros(nv, nfe)
    for mc in 1:nv, i in 1:nfe
        v_phys_mc = v_scaled[mc, i] * vs[mc]
        FO_L_new[mc,i] = (v_phys_mc - vlb[mc]) * alpha_L_new[mc,i]
        FO_U_new[mc,i] = (v_phys_mc - vub[mc]) * alpha_U_new[mc,i]
    end
    seed[:FO_L] = FO_L_new
    seed[:FO_U] = FO_U_new
    # los FO dinamicos quedan en 0 porque sus alphas estan en 0 desde NB1.
    # Forzamos a 0 explicitamente por consistencia con las constraints FO3..FO_atpm_c.
    if haskey(seed, :FO_upt);   fill!(seed[:FO_upt],   0.0); end
    if haskey(seed, :FO_prod);  fill!(seed[:FO_prod],  0.0); end
    if haskey(seed, :FO_aa);    fill!(seed[:FO_aa],    0.0); end
    if haskey(seed, :FO_pair);  fill!(seed[:FO_pair],  0.0); end
    if haskey(seed, :FO_ea);    fill!(seed[:FO_ea],    0.0); end
    if haskey(seed, :FO_atpm);  fill!(seed[:FO_atpm],  0.0); end

    @info "[dual-seed] lambda_/alpha_L/alpha_U/FO_* reconstruidos." *
          " residuo stationarity pre-alpha: max=$(round(maximum(res_before), sigdigits=4))" *
          " post-LSQ+activeset: max=$(round(maximum(res_after), sigdigits=4))" *
          " FO_L range=[$(round(minimum(FO_L_new), sigdigits=4)), $(round(maximum(FO_L_new), sigdigits=4))]" *
          " FO_U range=[$(round(minimum(FO_U_new), sigdigits=4)), $(round(maximum(FO_U_new), sigdigits=4))]"
    return seed
end

if seed !== nothing
    reconstruct_duals_from_seed!(seed)
end

# ─────────────────────────────────────────────────────────────────────────
# Diagnostico de residuos del seed: cuanta infactibilidad inicial trae,
# desglosada por familia de constraint LINEAL (las NL no se evaluan aqui).
# ─────────────────────────────────────────────────────────────────────────
function diagnose_seed_residuals(seed::Dict)
    if !(haskey(seed, :v) && haskey(seed, :c) && haskey(seed, :cdot) && haskey(seed, :hv))
        @warn "[seed-diag] faltan claves para diagnostico"
        return
    end
    v_s    = seed[:v]
    c_s    = seed[:c]
    cdot_s = seed[:cdot]
    hv_s   = seed[:hv]

    radau_colmat = [0.19681547722366  -0.06553542585020   0.02377097434822;
                    0.39442431473909   0.29207341166523  -0.04154875212600;
                    0.37640306270047   0.51248582618842   0.11111111111111]

    # ---- coll_c_0: c[l,1,j] = c0[l]/cs[l] + hv[1]*Σ colmat[j,k]*cdot[l,1,k] ----
    coll0_max = 0.0
    coll0_loc = (0,0,0)
    for l in 1:nc, j in 1:ncp
        rhs = c0[l]/cs[l]
        for k in 1:ncp
            rhs += hv_s[1] * radau_colmat[j,k] * cdot_s[l,1,k]
        end
        r = c_s[l,1,j] - rhs
        if abs(r) > coll0_max
            coll0_max = abs(r); coll0_loc = (l,1,j)
        end
    end

    # ---- coll_c_n: c[l,i,j] = c[l,i-1,ncp] + hv[i]*Σ colmat[j,k]*cdot[l,i,k] ----
    coll_n_max = 0.0
    coll_n_loc = (0,0,0)
    for l in 1:nc, i in 2:nfe, j in 1:ncp
        rhs = c_s[l,i-1,ncp]
        for k in 1:ncp
            rhs += hv_s[i] * radau_colmat[j,k] * cdot_s[l,i,k]
        end
        r = c_s[l,i,j] - rhs
        if abs(r) > coll_n_max
            coll_n_max = abs(r); coll_n_loc = (l,i,j)
        end
    end

    # ---- Sc balance: Σ S[mc,rx] * v[rx,i] * vs[rx] = 0 ----
    Sc_max = 0.0
    Sc_loc = (0,0)
    for i in 1:nfe
        for mc in 1:nm
            s = 0.0
            for rx in 1:nv
                s += S[mc,rx] * v_s[rx,i] * vs[rx]
            end
            if abs(s) > Sc_max
                Sc_max = abs(s); Sc_loc = (mc,i)
            end
        end
    end

    # ---- v_UB / v_LB slacks (negativos = infactible) ----
    vub_viol_max = 0.0
    vlb_viol_max = 0.0
    for rx in 1:nv, i in 1:nfe
        v_phys = v_s[rx,i] * vs[rx]
        slack_ub = vub[rx] - v_phys
        slack_lb = v_phys - vlb[rx]
        if slack_ub < -vub_viol_max; vub_viol_max = -slack_ub; end
        if slack_lb < -vlb_viol_max; vlb_viol_max = -slack_lb; end
    end

    # ---- MFE1 ----
    mfe1 = abs(sum(hv_s) - th)

    # ---- Stationarity Lagrangian residual con (lambda_, alpha) reconstruidos ----
    stat_max = 0.0
    if haskey(seed, :lambda_) && haskey(seed, :alpha_L) && haskey(seed, :alpha_U)
        lam_s = seed[:lambda_]
        aL_s  = seed[:alpha_L]
        aU_s  = seed[:alpha_U]
        # phase para d_phase
        phase_g_v = zeros(nfe)
        for i in 1:nfe
            cNfree_val = c_s[IDX_NFREE, i, ncp] * cs[IDX_NFREE]
            cNaa_val = 0.0
            for a in 1:length(AA_STATE_IDXS)
                cNaa_val += c_s[AA_STATE_IDXS[a], i, ncp] * cs[AA_STATE_IDXS[a]]
            end
            cNaa_val *= MW_N
            cNtot_val = cNfree_val + cNaa_val
            phase_g_v[i] = 1.0 / (1.0 + exp(-(cNtot_val - N_TOTAL_DEPLETION_THRESHOLD) / PHASE_SMOOTH_EPS))
        end
        for mc in 1:nv, i in 1:nfe
            d_mc = phase_g_v[i]*D_GROWTH[mc] + (1 - phase_g_v[i])*D_TURNOVER[mc] + phase_g_v[i]*D_AAUP[mc]
            v_phys = v_s[mc,i] * vs[mc]
            ST_lam = 0.0
            for r in 1:nm
                ST_lam += S[r,mc] * lam_s[r,i]
            end
            r_stat = d_mc + Q_REG*v_phys + aL_s[mc,i] + aU_s[mc,i] + ST_lam
            if abs(r_stat) > stat_max; stat_max = abs(r_stat); end
        end
    end

    println("\n[seed-diag] === Residuos LINEALES en el punto de partida ===")
    println("[seed-diag] coll_c_0   max |resid| = $(round(coll0_max,sigdigits=4))  loc=$coll0_loc")
    println("[seed-diag] coll_c_n   max |resid| = $(round(coll_n_max,sigdigits=4)) loc=$coll_n_loc")
    println("[seed-diag] Sc balance max |resid| = $(round(Sc_max,sigdigits=4))    loc=$Sc_loc")
    println("[seed-diag] v_UB violation max     = $(round(vub_viol_max,sigdigits=4))")
    println("[seed-diag] v_LB violation max     = $(round(vlb_viol_max,sigdigits=4))")
    println("[seed-diag] MFE1 (sum hv − th)     = $(round(mfe1,sigdigits=4))")
    println("[seed-diag] Stationarity (LSQ)     = $(round(stat_max,sigdigits=4))")
    println("[seed-diag] FO_L max |val|         = $(round(maximum(abs.(seed[:FO_L])),sigdigits=4))")
    println("[seed-diag] FO_U max |val|         = $(round(maximum(abs.(seed[:FO_U])),sigdigits=4))")
    println()
end

if seed !== nothing
    diagnose_seed_residuals(seed)
end


include("pFBA_KKT_flux_Zenteno_vargam_simultaneous_v2.jl")

@assert @isdefined(c0) "c0 no definido. Ejecutar celda de Malla/estados/c0."

println("Lanzando MPCC v2 (Scholtes homotopy)...")
println("  nc=$nc, nfe=$nfe, ncp=$ncp, th=$th")
println("  seed = ", seed === nothing ? "NO" : "SI ($(length(keys(seed))) claves)")
println("  eps_schedule = ", EPS_SCHEDULE)

result = pFBA_KKT_flux_Zenteno_vargam_simultaneous_v2(
    c0;
    eps_flux = EPS_FLUX,
    apply_product_caps = APPLY_PRODUCT_CAPS,
    seed = seed,
    eps_schedule = EPS_SCHEDULE,
    verbose_homotopy = true,
)

cStar, vStar, cdStar, lStar, alLStar, alUStar,
    auStar, apStar, aaaStar, aprStar, aeaStar, aatpmStar, hStar, diagStar = result

println("\nResolucion terminada.")
println("Status = ", diagStar.solver)


# Reconstruir grilla temporal y trayectoria de estados (en unidades FISICAS)
radau_roots = [0.15505102572168, 0.64494897427832, 1.0]
t_fe = [sum(hStar[1:i-1]) for i in 1:nfe]
t_col = Float64[]
for i in 1:nfe, j in 1:ncp
    push!(t_col, t_fe[i] + radau_roots[j] * hStar[i])
end

println("h*: min=$(minimum(hStar))  max=$(maximum(hStar))  sum=$(sum(hStar))")
println("nfe*ncp puntos temporales = ", length(t_col))

# Trayectoria en unidades fisicas: el builder YA desescala cStar (c_phys = c/cs * cs).
# No volver a multiplicar aqui (bug previo de doble-escalamiento).
Y = zeros(nc, length(t_col))
for i in 1:nfe, j in 1:ncp
    for s in 1:nc
        Y[s, (i-1)*ncp + j] = cStar[s, i, j]
    end
end
println("Trayectoria shape: ", size(Y))


plot_rows = [
    (IDX_X,     "X (gDW/L)"),
    (IDX_NFREE, "N libre (gN/L)"),
    (IDX_G,     "Glucosa (g/L)"),
    (IDX_F,     "Fructosa (g/L)"),
    (IDX_E,     "Etanol (g/L)"),
    (IDX_PROT,  "Prot (g/L)"),
    (IDX_CARB,  "Carb (g/L)"),
    (IDX_NREC,  "N_rec (gN/L)"),
]

plts = [plot(t_col, Y[s,:], title=ttl, lw=2, xlabel="t [h]", legend=false)
        for (s, ttl) in plot_rows]
plot(plts..., layout=(length(plts), 1), size=(700, 200*length(plts)))

# Trayectoria de estados (unidades fisicas)
open(joinpath(OUT_DIR, "xk_simultaneous_v2.csv"), "w") do f
    for i in 1:size(Y,2)
        println(f, join(Y[:,i], ","))
    end
end

# Flujos: el builder YA desescala vStar. No re-multiplicar.
open(joinpath(OUT_DIR, "v_simultaneous_v2.csv"), "w") do f
    for i in 1:nfe
        v_phys = [vStar[rx,i] for rx in 1:nv]
        println(f, join(v_phys, ","))
    end
end

# Tiempos
open(joinpath(OUT_DIR, "t_simultaneous_v2.csv"), "w") do f
    for t in t_col; println(f, t); end
end

println("CSV v2 exportados a ", OUT_DIR)