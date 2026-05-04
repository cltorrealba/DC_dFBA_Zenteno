# ======================================================================================
# pFBA_KKT_flux_Zenteno_vargam_simultaneous.jl
# --------------------------------------------------------------------------------------
# Simultaneous direct-collocation / MPCC update aligned with Notebook 1's new formulation,
# but with the key IPOPT-stabilizing approximations described below:
#
# 1) Stoichiometric matrix S is kept FIXED.
# 2) Variable GAM is injected as an exogenous ATP burden schedule GAM_EXTRA_FE[i],
#    typically exported from Notebook 1, instead of modifying biomass stoichiometry
#    inside the NLP.
# 3) The phase switch BIOMASS -> TURNOVER_ATPM is represented smoothly via N_total.
# 4) The exact two-stage pFBA secondary problem is relaxed to:
#       - dynamic LP objective coefficients (biomass vs ATPM)
#       - optional AA uptake incentive during growth
#       - tiny Tikhonov regularization on fluxes
# 5) Pairwise acetate-ester couplings and ethyl-acetate soft coupling remain explicit.
#
# Expected globals to be defined by the notebook:
#   S, vlb, vub, RXN_IDS, MET_IDS, nm, nv
#   nfe, ncp, nc, th, hm, var_h
#   vs, cs
#   obj, glu, fru, eth, o2, IDX_ATPM, IDX_PROT_RXN
#   IDX_X, IDX_NFREE, IDX_G, IDX_F, IDX_E, IDX_O2, IDX_PROT, IDX_CARB
#   KINETIC_N_SOURCE_IDXS, UPTAKE_IDXS, PRODUCT_IDXS
#   N_atoms_vec, N_frac, SELECT_UPTAKE, SELECT_PRODUCT
#   AA_STATE_IDXS, AA_UPTAKE_IDXS, AA_ALPHA_VEC, SELECT_AAUPTAKE
#   AROMA_STATE_IDXS, AROMA_RXN_IDXS, AROMA_MW_VEC
#   PAIRWISE_ESTER_IDXS, PAIRWISE_ALCOHOL_IDXS, PAIRWISE_PHI
#   EA_SOFT_ESTER_IDX, EA_SOFT_ALCOHOL_IDX, PHI_ETHYL_ACETATE_STATIC
#   GAM_EXTRA_FE
#   kinetic constants and helper functions used below.
# ======================================================================================

using JuMP
using Ipopt
using LinearAlgebra
using Statistics
using SparseArrays
using Libdl
import MathOptInterface as MOI

_softplus(x, ϵ) = 0.5 * (x + sqrt(x * x + ϵ * ϵ))
_sigmoid(x) = 1.0 / (1.0 + exp(-x))

const HSL_LINEAR_SOLVERS = Set(["ma27", "ma57", "ma77", "ma86", "ma97"])

const HSL_JLL_AVAILABLE = let
    try
        import HSL_jll
        true
    catch
        false
    end
end

function _hsl_lib_path_from_env()
    for key in ("IPOPT_HSLLIB", "HSL_DLL_PATH", "COINHSL_DLL_PATH")
        candidate = strip(get(ENV, key, ""))
        if !isempty(candidate) && isfile(candidate)
            return String(candidate), key
        end
    end
    return nothing, nothing
end

function _is_hsl_library_functional(lib_path::AbstractString)
    handle = nothing
    try
        handle = Libdl.dlopen(lib_path)
        try
            sym = Libdl.dlsym(handle, :LIBHSL_isfunctional)
            return ccall(sym, Cint, ()) != 0
        catch
            @warn "La librería HSL no exporta LIBHSL_isfunctional; se asume funcional si carga correctamente." lib_path
            return true
        end
    catch err
        @warn "No se pudo cargar la librería HSL indicada." lib_path exception = (err, catch_backtrace())
        return false
    finally
        handle !== nothing && try Libdl.dlclose(handle) catch end
    end
end

function _configure_ipopt_linear_solver(requested_solver::AbstractString)
    solver = lowercase(String(requested_solver))
    hsl_lib_path = nothing
    hsl_functional = false

    if solver in HSL_LINEAR_SOLVERS
        env_hsl_path, env_key = _hsl_lib_path_from_env()
        if env_hsl_path !== nothing
            env_hsl_ok = _is_hsl_library_functional(env_hsl_path)
            if env_hsl_ok
                hsl_lib_path = env_hsl_path
                hsl_functional = true
                println("[ipopt] hsllib desde ", env_key, " = ", env_hsl_path)
            else
                @warn "La librería HSL indicada en entorno no es funcional." env_key env_hsl_path
            end
        end

        if !hsl_functional && !HSL_JLL_AVAILABLE
            @warn "No se pudo cargar HSL_jll y no hay hsllib funcional en entorno; se usará MUMPS." requested_solver = solver
            solver = "mumps"
        elseif !hsl_functional && HSL_JLL_AVAILABLE
            try
            hsl_functional = @ccall HSL_jll.libhsl.LIBHSL_isfunctional()::Bool
            if hsl_functional
                hsl_lib_path = String(HSL_jll.libhsl_path)
            else
                @warn "HSL_jll está instalado pero la libHSL cargada no es funcional; se usará MUMPS." requested_solver = solver
                solver = "mumps"
            end
            catch err
                @warn "No se pudo inicializar HSL_jll; se usará MUMPS." requested_solver = solver exception = (err, catch_backtrace())
                solver = "mumps"
            end
        end

        if solver == "ma86" && hsl_functional
            try
                blas_threads = LinearAlgebra.BLAS.get_num_threads()
                if blas_threads != 1
                    LinearAlgebra.BLAS.set_num_threads(1)
                    println("[ipopt] BLAS threads ajustados a 1 para evitar oversubscription con ma86.")
                end
            catch err
                @warn "No fue posible ajustar los threads de BLAS para ma86." exception = (err, catch_backtrace())
            end
        end
    end

    return solver, hsl_lib_path, hsl_functional
end

function _build_ipopt_attributes()
    requested_solver = lowercase(String(IPOPT_LINEAR_SOLVER))
    active_solver, hsl_lib_path, hsl_functional = _configure_ipopt_linear_solver(requested_solver)
    use_dual_warm_start = lowercase(get(ENV, "IPOPT_USE_DUAL_WARM_START", "false")) in ("1", "true", "yes", "on")
    warm_start_init_point = lowercase(get(ENV, "IPOPT_WARM_START_INIT_POINT", "yes")) in ("1", "true", "yes", "on") ? "yes" : "no"

    attrs = Pair{String,Any}[
        "linear_solver" => active_solver,
        "print_level" => IPOPT_PRINT_LEVEL,
        "tol" => IPOPT_TOL,
        "acceptable_tol" => IPOPT_ACCEPTABLE_TOL,
        "acceptable_iter" => IPOPT_ACCEPTABLE_ITER,
        "max_iter" => IPOPT_MAX_ITER,
        "mu_strategy" => "adaptive",
        "mu_init" => 1e-2,
        "mu_min" => 1e-8,
        "constr_viol_tol" => IPOPT_CONSTR_VIOL_TOL,
        "compl_inf_tol" => IPOPT_COMPL_INF_TOL,
        "bound_relax_factor" => 1e-8,
        "warm_start_init_point" => warm_start_init_point,
        "mumps_mem_percent" => try parse(Int, get(ENV, "IPOPT_MUMPS_MEM_PERCENT", "20")) catch; 20 end,
    ]

    if hsl_lib_path !== nothing
        push!(attrs, "hsllib" => hsl_lib_path)
    end

    println("[ipopt] linear_solver requested=", requested_solver, " active=", active_solver, " hsl_ready=", hsl_functional)
    println("[ipopt] dual_warm_start=", use_dual_warm_start, " | warm_start_init_point=", warm_start_init_point)
    hsl_lib_path !== nothing && println("[ipopt] hsllib = ", hsl_lib_path)

    return attrs
end

function pFBA_KKT_flux_Zenteno_vargam_simultaneous(
    c0;
    eps_flux::Float64 = 0.0,
    apply_product_caps::Bool = false,
    # ── primal ──────────────────────────────────────────────────────────────
    warm_c::Union{Nothing, Array{Float64,3}} = nothing,
    warm_v::Union{Nothing, Matrix{Float64}} = nothing,
    warm_h::Union{Nothing, Vector{Float64}} = nothing,
    warm_cdot::Union{Nothing, Array{Float64,3}} = nothing,
    # ── dual / Lagrange multipliers ─────────────────────────────────────────
    warm_lambda::Union{Nothing, Matrix{Float64}} = nothing,
    warm_alL::Union{Nothing, Matrix{Float64}} = nothing,
    warm_alU::Union{Nothing, Matrix{Float64}} = nothing,
    warm_alupt::Union{Nothing, Matrix{Float64}} = nothing,
    warm_alprod::Union{Nothing, Matrix{Float64}} = nothing,
    warm_alaa::Union{Nothing, Matrix{Float64}} = nothing,
    warm_alpair::Union{Nothing, Matrix{Float64}} = nothing,
    warm_alea::Union{Nothing, Vector{Float64}} = nothing,
    warm_alatpm::Union{Nothing, Vector{Float64}} = nothing,
)
    @assert length(c0) == nc "c0 debe tener longitud nc."
    @assert nc >= 8 "Este modelo simultáneo asume al menos 8 estados base."

    # ---------- metadata sizes ----------
    N_UP   = length(UPTAKE_IDXS)
    N_PROD = length(PRODUCT_IDXS)
    N_AA   = length(AA_UPTAKE_IDXS)
    N_ARO  = length(AROMA_RXN_IDXS)
    N_PAIR = length(PAIRWISE_ESTER_IDXS)
    HAS_EA_SOFT = (EA_SOFT_ESTER_IDX > 0) && (EA_SOFT_ALCOHOL_IDX > 0) && (PHI_ETHYL_ACETATE_STATIC > 0.0)

    # ---------- sparse structural maps (critical for IPOPT memory) ----------
    Ssp = sparse(S)
    S_I, S_J, S_V = findnz(Ssp)

    S_ROW_COLS = [Int[] for _ in 1:nm]
    S_ROW_VALS = [Float64[] for _ in 1:nm]
    S_COL_ROWS = [Int[] for _ in 1:nv]
    S_COL_VALS = [Float64[] for _ in 1:nv]

    for q in eachindex(S_V)
        r = S_I[q]
        c = S_J[q]
        vq = Float64(S_V[q])
        push!(S_ROW_COLS[r], c)
        push!(S_ROW_VALS[r], vq)
        push!(S_COL_ROWS[c], r)
        push!(S_COL_VALS[c], vq)
    end

    UPTAKE_ACTIVE = [Int[] for _ in 1:nv]
    for k in 1:N_UP
        mc = UPTAKE_IDXS[k]
        (1 <= mc <= nv) && push!(UPTAKE_ACTIVE[mc], k)
    end

    PRODUCT_ACTIVE = [Int[] for _ in 1:nv]
    for k in 1:N_PROD
        mc = PRODUCT_IDXS[k]
        (1 <= mc <= nv) && push!(PRODUCT_ACTIVE[mc], k)
    end

    AA_ACTIVE = [Int[] for _ in 1:nv]
    for a in 1:N_AA
        mc = AA_UPTAKE_IDXS[a]
        (1 <= mc <= nv) && push!(AA_ACTIVE[mc], a)
    end

    PAIR_ACTIVE_IDXS = [Int[] for _ in 1:nv]
    PAIR_ACTIVE_COEFFS = [Float64[] for _ in 1:nv]
    for p in 1:N_PAIR
        mc_alc = PAIRWISE_ALCOHOL_IDXS[p]
        mc_est = PAIRWISE_ESTER_IDXS[p]
        if 1 <= mc_alc <= nv
            push!(PAIR_ACTIVE_IDXS[mc_alc], p)
            push!(PAIR_ACTIVE_COEFFS[mc_alc], Float64(PAIRWISE_PHI[p]))
        end
        if 1 <= mc_est <= nv
            push!(PAIR_ACTIVE_IDXS[mc_est], p)
            push!(PAIR_ACTIVE_COEFFS[mc_est], -1.0)
        end
    end

    EA_STATIONARITY_COEFF = zeros(Float64, nv)
    if HAS_EA_SOFT
        if 1 <= EA_SOFT_ALCOHOL_IDX <= nv
            EA_STATIONARITY_COEFF[EA_SOFT_ALCOHOL_IDX] += PHI_ETHYL_ACETATE_STATIC
        end
        if 1 <= EA_SOFT_ESTER_IDX <= nv
            EA_STATIONARITY_COEFF[EA_SOFT_ESTER_IDX] -= 1.0
        end
    end

    println("[sparse] nnz(S) = ", length(S_V), " | density = ", round(length(S_V) / (nm * nv), digits=6))

    # ---------- Radau IIA ----------
    colmat = [0.19681547722366  -0.06553542585020   0.02377097434822;
              0.39442431473909   0.29207341166523  -0.04154875212600;
              0.37640306270047   0.51248582618842   0.11111111111111]
    radau  = [0.15505102572168, 0.64494897427832, 1.0]

    # ---------- Model ----------
    ipopt_attrs = _build_ipopt_attributes()
    m = Model(optimizer_with_attributes(Ipopt.Optimizer, ipopt_attrs...))

    # ---------- Variables ----------
    @variables(m, begin
        c[1:nc, 1:nfe, 1:ncp] >= 0.0
        cdot[1:nc, 1:nfe, 1:ncp]
        v[1:nv, 1:nfe]

        lambda_[1:nm, 1:nfe]

        alpha_U[1:nv, 1:nfe] >= 0.0
        alpha_L[1:nv, 1:nfe] <= 0.0

        alpha_upt[1:N_UP, 1:nfe] <= 0.0
        alpha_prod[1:N_PROD, 1:nfe] >= 0.0
        alpha_aa[1:N_AA, 1:nfe] <= 0.0
        alpha_pair[1:N_PAIR, 1:nfe] >= 0.0
        alpha_ea[1:nfe] >= 0.0
        alpha_atpm_floor[1:nfe] >= 0.0

        FO_L[1:nv, 1:nfe]
        FO_U[1:nv, 1:nfe]
        FO_upt[1:N_UP, 1:nfe]
        FO_prod[1:N_PROD, 1:nfe]
        FO_aa[1:N_AA, 1:nfe]
        FO_pair[1:N_PAIR, 1:nfe]
        FO_ea[1:nfe]
        FO_atpm[1:nfe]

        hv[1:nfe] >= 0.0
    end)

    # ---------- starts ----------
    # IMPORTANTE: warm_c y warm_v deben estar en ESCALA INTERNA del modelo (divididas por cs/vs)
    _use_wc = (warm_c !== nothing) && size(warm_c) == (nc, nfe, ncp)
    _use_wv = (warm_v !== nothing) && size(warm_v) == (nv, nfe)
    _use_wh = (warm_h !== nothing) && length(warm_h) == nfe
    if any((_use_wc, _use_wv, _use_wh))
        println("[WS] Warm start aplicado: c=$(_use_wc), v=$(_use_wv), h=$(_use_wh)")
        if _use_wc
            c_sample = warm_c[1,1,1]
            c_ref = c0[1] / cs[1]
            c_ratio = c_sample > 0 ? c_sample / max(c_ref, 1e-9) : 1.0
            if c_ratio > 2.0 || c_ratio < 0.5
                @warn "[WS] Escala de warm_c parece incorrecta (ratio vs default=$c_ratio). Asegurate de que warm_c esté dividido por cs."
            end
        end
        if _use_wv
            v_sample = warm_v[max(1, nv÷2), 1]
            v_ref = 0.0
            if abs(v_sample) > max(maximum(abs.(vub)), maximum(abs.(vlb))) * 0.1
                @warn "[WS] Escala de warm_v parece anómala. Asegurate de que warm_v esté dividido por vs."
            end
        end
    end

    for i in 1:nfe
        set_start_value(hv[i], _use_wh ? warm_h[i] : hm[i])
        for j in 1:ncp
            for s in 1:nc
                set_start_value(c[s, i, j], _use_wc ? warm_c[s, i, j] : c0[s] / cs[s])
            end
        end
        for rx in 1:nv
            set_start_value(v[rx, i], _use_wv ? warm_v[rx, i] : 0.0)
        end
    end

    # ----- derivative warm start (cdot) -----
    # IMPORTANTE: warm_cdot debe estar en ESCALA INTERNA (d/dt de c escalado, SIN dividir por cs)
    if warm_cdot !== nothing && size(warm_cdot) == (nc, nfe, ncp)
        for i in 1:nfe, j in 1:ncp, s in 1:nc
            set_start_value(cdot[s, i, j], warm_cdot[s, i, j])
        end
    end

    # ----- dual / multiplier warm starts -----
    use_dual_warm_start = lowercase(get(ENV, "IPOPT_USE_DUAL_WARM_START", "false")) in ("1", "true", "yes", "on")
    if use_dual_warm_start
        if warm_lambda !== nothing && size(warm_lambda) == (nm, nfe)
            for i in 1:nfe, r in 1:nm
                set_start_value(lambda_[r, i], warm_lambda[r, i])
            end
        end
        if warm_alL !== nothing && size(warm_alL) == (nv, nfe)
            for i in 1:nfe, mc in 1:nv
                set_start_value(alpha_L[mc, i], min(0.0, warm_alL[mc, i]))
            end
        end
        if warm_alU !== nothing && size(warm_alU) == (nv, nfe)
            for i in 1:nfe, mc in 1:nv
                set_start_value(alpha_U[mc, i], max(0.0, warm_alU[mc, i]))
            end
        end
        if N_UP > 0 && warm_alupt !== nothing && size(warm_alupt) == (N_UP, nfe)
            for i in 1:nfe, k in 1:N_UP
                set_start_value(alpha_upt[k, i], min(0.0, warm_alupt[k, i]))
            end
        end
        if N_PROD > 0 && warm_alprod !== nothing && size(warm_alprod) == (N_PROD, nfe)
            for i in 1:nfe, k in 1:N_PROD
                set_start_value(alpha_prod[k, i], max(0.0, warm_alprod[k, i]))
            end
        end
        if N_AA > 0 && warm_alaa !== nothing && size(warm_alaa) == (N_AA, nfe)
            for i in 1:nfe, a in 1:N_AA
                set_start_value(alpha_aa[a, i], min(0.0, warm_alaa[a, i]))
            end
        end
        if N_PAIR > 0 && warm_alpair !== nothing && size(warm_alpair) == (N_PAIR, nfe)
            for i in 1:nfe, p in 1:N_PAIR
                set_start_value(alpha_pair[p, i], max(0.0, warm_alpair[p, i]))
            end
        end
        if warm_alea !== nothing && length(warm_alea) == nfe
            for i in 1:nfe
                set_start_value(alpha_ea[i], max(0.0, warm_alea[i]))
            end
        end
        if warm_alatpm !== nothing && length(warm_alatpm) == nfe
            for i in 1:nfe
                set_start_value(alpha_atpm_floor[i], max(0.0, warm_alatpm[i]))
            end
        end
        let _ws_dual = String[]
            warm_cdot   !== nothing && push!(_ws_dual, "cdot")
            warm_lambda !== nothing && push!(_ws_dual, "\u03bb")
            warm_alL    !== nothing && push!(_ws_dual, "\u03b1L")
            warm_alU    !== nothing && push!(_ws_dual, "\u03b1U")
            warm_alupt  !== nothing && push!(_ws_dual, "\u03b1UPT")
            warm_alprod !== nothing && push!(_ws_dual, "\u03b1PROD")
            warm_alaa   !== nothing && push!(_ws_dual, "\u03b1AA")
            warm_alpair !== nothing && push!(_ws_dual, "\u03b1PAIR")
            warm_alea   !== nothing && push!(_ws_dual, "\u03b1EA")
            warm_alatpm !== nothing && push!(_ws_dual, "\u03b1ATPM")
            !isempty(_ws_dual) && println("[WS] Dual/derivadas: ", join(_ws_dual, ", "))
        end
    else
        println("[WS] Dual warm-start desactivado (IPOPT_USE_DUAL_WARM_START=false).")
    end

    c0s = [c0[i] / cs[i] for i in 1:nc]

    # ---------- register helper functions ----------
    JuMP.register(m, :dynamic_temperature, 1, dynamic_temperature; autodiff = true)
    JuMP.register(m, :death_rate_T, 2, death_rate_T; autodiff = true)
    JuMP.register(m, :smooth_injection, 4, smooth_injection; autodiff = true)
    JuMP.register(m, :softplus_jl, 2, _softplus; autodiff = true)
    JuMP.register(m, :sigmoid_jl, 1, _sigmoid; autodiff = true)

    # ---------- time map ----------
    @NLexpression(m, t_fe[i=1:nfe], sum(hv[k] for k in 1:i-1))
    @NLexpression(m, t_loc[i=1:nfe, j=1:ncp], t_fe[i] + radau[j] * hv[i])

    # ---------- scaled-to-physical state aliases ----------
    @NLexpressions(m, begin
        cX[i=1:nfe, j=1:ncp], c[IDX_X, i, j] * cs[IDX_X]
        cNf[i=1:nfe, j=1:ncp], c[IDX_NFREE, i, j] * cs[IDX_NFREE]
        cG[i=1:nfe, j=1:ncp], c[IDX_G, i, j] * cs[IDX_G]
        cF[i=1:nfe, j=1:ncp], c[IDX_F, i, j] * cs[IDX_F]
        cE[i=1:nfe, j=1:ncp], c[IDX_E, i, j] * cs[IDX_E]
        cO2[i=1:nfe, j=1:ncp], c[IDX_O2, i, j] * cs[IDX_O2]
        cProt[i=1:nfe, j=1:ncp], c[IDX_PROT, i, j] * cs[IDX_PROT]
        cCarb[i=1:nfe, j=1:ncp], c[IDX_CARB, i, j] * cs[IDX_CARB]
        cNaa[i=1:nfe, j=1:ncp], MW_N * sum(c[AA_STATE_IDXS[a], i, j] * cs[AA_STATE_IDXS[a]] for a in 1:N_AA)
        cNtot[i=1:nfe, j=1:ncp], cNf[i,j] + cNaa[i,j]
        phase_g[i=1:nfe, j=1:ncp], sigmoid_jl((cNtot[i,j] - N_TOTAL_DEPLETION_THRESHOLD) / PHASE_SMOOTH_EPS)
        phase_t[i=1:nfe, j=1:ncp], 1.0 - phase_g[i,j]
    end)

    @NLexpression(m, phase_g_fe[i=1:nfe], phase_g[i, ncp])
    @NLexpression(m, phase_t_fe[i=1:nfe], 1.0 - phase_g_fe[i])

    # ---------- temperature / kinetic envelopes ----------
    @NLexpression(m, T_loc[i=1:nfe, j=1:ncp], dynamic_temperature(t_loc[i,j]))
    @NLexpression(m, mu_T_ij[i=1:nfe, j=1:ncp],
        exp(59453.0 * (T_loc[i,j] - 300.0) / (300.0 * R * T_loc[i,j])))
    @NLexpression(m, Kg_T_ij[i=1:nfe, j=1:ncp],
        exp(46055.0 * (T_loc[i,j] - 293.15) / (293.15 * R * T_loc[i,j])))
    @NLexpression(m, b_T_ij[i=1:nfe, j=1:ncp],
        exp(11000.0 * (T_loc[i,j] - 296.15) / (296.15 * R * T_loc[i,j])))
    @NLexpression(m, mrate_ij[i=1:nfe, j=1:ncp],
        MRATE_0 * exp(37681.0 * (T_loc[i,j] - 293.30) / (293.30 * R * T_loc[i,j])))

    @NLexpression(m, mu_j[i=1:nfe, j=1:ncp],
        MU0 * mu_T_ij[i,j] * (cNf[i,j] / (cNf[i,j] + Kn0_nom * Kg_T_ij[i,j] + EPS)))
    @NLexpression(m, betaG_j[i=1:nfe, j=1:ncp],
        betaG0_nom * b_T_ij[i,j] *
        (cG[i,j] / (cG[i,j] + Kg0_nom * Kg_T_ij[i,j] + EPS)) *
        (Kie0_nom * Kg_T_ij[i,j] / (cE[i,j] + Kie0_nom * Kg_T_ij[i,j] + EPS)))
    @NLexpression(m, betaF_j[i=1:nfe, j=1:ncp],
        betaF0_nom * b_T_ij[i,j] *
        (cF[i,j] / (cF[i,j] + Kf0_nom * Kg_T_ij[i,j] + EPS)) *
        (Kig0_nom * Kg_T_ij[i,j] / (cG[i,j] + Kig0_nom * Kg_T_ij[i,j] + EPS)) *
        (Kie0_nom * Kg_T_ij[i,j] / (cE[i,j] + Kie0_nom * Kg_T_ij[i,j] + EPS)))
    @NLexpression(m, Kd_j[i=1:nfe, j=1:ncp], death_rate_T(cE[i,j], T_loc[i,j]))

    @NLexpression(m, injection_rate[i=1:nfe, j=1:ncp],
        smooth_injection(t_loc[i,j], T_INJ_1, DOSE_1, WIDTH_1) +
        smooth_injection(t_loc[i,j], T_INJ_2, DOSE_2, WIDTH_2))

    # ---------- kinetic macro-limits (same spirit as old model) ----------
    @NLexpressions(m, begin
        vx[i=1:nfe, j=1:ncp], mu_j[i,j]
        vg[i=1:nfe, j=1:ncp], mu_j[i,j] / YXG_nom +
                               betaG_j[i,j] / YEG +
                               mrate_ij[i,j] * (cG[i,j] / (cG[i,j] + cF[i,j] + EPS))
        vf[i=1:nfe, j=1:ncp], mu_j[i,j] / YXF_nom +
                               betaF_j[i,j] / YEF +
                               mrate_ij[i,j] * (cF[i,j] / (cG[i,j] + cF[i,j] + EPS))
        vn[i=1:nfe, j=1:ncp], mu_j[i,j] / YXN
        ve[i=1:nfe, j=1:ncp], (betaG_j[i,j] + betaF_j[i,j]) / MW_ETH

        L_uptake[k=1:N_UP, i=1:nfe, j=1:ncp],
            IS_GLU[k] * (vg[i,j] / MW_GLU) +
            IS_FRU[k] * (vf[i,j] / MW_FRU) +
            IS_NIT[k] * vn[i,j] * N_frac[k]

        L_product[k=1:N_PROD, i=1:nfe, j=1:ncp],
            IS_ETH_prod[k] * ve[i,j] + IS_OBJ_prod[k] * vx[i,j]
    end)

    PROD_CAP_ON = apply_product_caps ? 1.0 : 0.0
    PROD_CAP_BIG = 1.0e6
    @NLexpression(m, L_product_cap[k=1:N_PROD, i=1:nfe],
        L_product[k,i,ncp] + (1.0 - PROD_CAP_ON) * PROD_CAP_BIG)

    # ---------- AA uptake / turnover caps ----------
    if N_AA > 0
        @NLexpression(m, Q_AA_growth[a=1:N_AA, i=1:nfe],
            K_AA_UPTAKE_GROWTH * (c[AA_STATE_IDXS[a], i, ncp] * cs[AA_STATE_IDXS[a]]) / (cX[i,ncp] + EPS))
        @NLexpression(m, Q_AA_turn[a=1:N_AA, i=1:nfe],
            TURNOVER_LAMBDA * cProt[i,ncp] * AA_ALPHA_VEC[a] / (XA_FRACTION * cX[i,ncp] + EPS))
        @NLexpression(m, Q_AA_cap[a=1:N_AA, i=1:nfe],
            phase_g_fe[i] * Q_AA_growth[a,i] + phase_t_fe[i] * Q_AA_turn[a,i])
    end

    # ---------- ATP maintenance / variable GAM burden ----------
    @NLexpression(m, ATPM_floor_rhs[i=1:nfe],
        phase_t_fe[i] * ATPM_LB_NO_GROWTH + GAM_EXTRA_FE[i] * softplus_jl(v[obj,i] * vs[obj], SOFTPLUS_V_EPS))

    # ---------- Dynamic objective coefficients ----------
    # Growth: maximize biomass. Turnover: maximize ATPM. Growth also mildly rewards AA uptake.
    @NLexpression(m, d_phase[mc=1:nv, i=1:nfe],
        phase_g_fe[i] * D_GROWTH[mc] + phase_t_fe[i] * D_TURNOVER[mc] + phase_g_fe[i] * D_AAUP[mc])

    # ---------- Outer objective ----------
    @NLobjective(m, Min,
        sum(
            sum(-PHI_L * FO_L[mc,i] - PHI_U * FO_U[mc,i] for mc in 1:nv) +
            sum(PHI_UPT * FO_upt[k,i] for k in 1:N_UP) +
            sum(-PHI_PROD * FO_prod[k,i] for k in 1:N_PROD) +
            sum(PHI_AA * FO_aa[a,i] for a in 1:N_AA) +
            sum(-PHI_PAIR * FO_pair[p,i] for p in 1:N_PAIR) +
            (-PHI_EA * FO_ea[i]) +
            (-PHI_ATPM * FO_atpm[i]) +
            FLUX_SMOOTH_WEIGHT * sum((v[mc,i] * vs[mc])^2 for mc in 1:nv)
        for i in 1:nfe)
    )

    # ---------- linear constraints ----------
    @constraints(m, begin
        coll_c_0[l=1:nc, j=1:ncp],
            c[l,1,j] == c0s[l] + hv[1] * sum(colmat[j,k] * cdot[l,1,k] for k in 1:ncp)

        coll_c_n[l=1:nc, i=2:nfe, j=1:ncp],
            c[l,i,j] == c[l,i-1,ncp] + hv[i] * sum(colmat[j,k] * cdot[l,i,k] for k in 1:ncp)

        v_UB[mc=1:nv, i=1:nfe],  v[mc,i] * vs[mc] - vub[mc] <= 0
        v_LB[mc=1:nv, i=1:nfe], -v[mc,i] * vs[mc] + vlb[mc] <= 0

        MFE1, sum(hv[i] for i in 1:nfe) == th
        MFE3[i=1:nfe], hv[i] >= (1.0 - var_h) * hm[1]
        MFE4[i=1:nfe], hv[i] <= (1.0 + var_h) * hm[1]

        flux_smooth_pos[mc=1:nv, i=2:nfe],  v[mc,i] - v[mc,i-1] <= eps_flux
        flux_smooth_neg[mc=1:nv, i=2:nfe],  v[mc,i-1] - v[mc,i] <= eps_flux
    end)

    Sc = Matrix{ConstraintRef}(undef, nm, nfe)
    for mc in 1:nm
        cols = S_ROW_COLS[mc]
        vals = S_ROW_VALS[mc]
        for i in 1:nfe
            Sc[mc,i] = @constraint(m, sum(vals[t] * v[cols[t], i] * vs[cols[t]] for t in eachindex(cols)) == 0)
        end
    end

    # ---------- stationarity (w.r.t. fluxes) ----------
    lagr_stationarity = Matrix{ConstraintRef}(undef, nv, nfe)
    for mc in 1:nv
        up_terms   = UPTAKE_ACTIVE[mc]
        prod_terms = PRODUCT_ACTIVE[mc]
        aa_terms   = AA_ACTIVE[mc]
        pair_idxs  = PAIR_ACTIVE_IDXS[mc]
        pair_coeff = PAIR_ACTIVE_COEFFS[mc]
        srows      = S_COL_ROWS[mc]
        svals      = S_COL_VALS[mc]
        ea_coeff   = EA_STATIONARITY_COEFF[mc]
        obj_coeff  = (mc == obj) ? 1.0 : 0.0
        atpm_coeff = (mc == IDX_ATPM) ? -1.0 : 0.0

        for i in 1:nfe
            lagr_stationarity[mc,i] = @NLconstraint(m,
                d_phase[mc,i] +
                Q_REG * (v[mc,i] * vs[mc]) +
                alpha_L[mc,i] + alpha_U[mc,i] +
                sum(alpha_upt[up_terms[t],i] for t in 1:length(up_terms)) +
                sum(alpha_prod[prod_terms[t],i] for t in 1:length(prod_terms)) +
                sum(alpha_aa[aa_terms[t],i] for t in 1:length(aa_terms)) +
                sum(pair_coeff[t] * alpha_pair[pair_idxs[t],i] for t in 1:length(pair_idxs)) +
                ea_coeff * alpha_ea[i] +
                (obj_coeff * GAM_EXTRA_FE[i] + atpm_coeff) * alpha_atpm_floor[i] +
                sum(svals[t] * lambda_[srows[t],i] for t in 1:length(srows)) == 0
            )
        end
    end

    # ---------- nonlinear dynamic + coupling constraints ----------
    @NLconstraints(m, begin
        # macro balances
        ode_X[i=1:nfe, j=1:ncp],
            cdot[IDX_X,i,j] == ((phase_g[i,j] * softplus_jl(v[obj,i] * vs[obj], SOFTPLUS_V_EPS) - Kd_j[i,j]) * cX[i,j]) / cs[IDX_X]

        ode_Nfree[i=1:nfe, j=1:ncp],
            cdot[IDX_NFREE,i,j] ==
                (- MW_N * phase_g[i,j] *
                    sum(IS_NIT[k] * N_atoms_vec[UPTAKE_IDXS[k]] *
                        (-v[UPTAKE_IDXS[k],i] * vs[UPTAKE_IDXS[k]]) for k in 1:N_UP) * cX[i,j]
                 + injection_rate[i,j]) / cs[IDX_NFREE]

        ode_G[i=1:nfe, j=1:ncp],
            cdot[IDX_G,i,j] == -MW_GLU * (-v[glu,i] * vs[glu]) * cX[i,j] / cs[IDX_G]

        ode_F[i=1:nfe, j=1:ncp],
            cdot[IDX_F,i,j] == -MW_FRU * (-v[fru,i] * vs[fru]) * cX[i,j] / cs[IDX_F]

        ode_E[i=1:nfe, j=1:ncp],
            cdot[IDX_E,i,j] == MW_ETH * softplus_jl(v[eth,i] * vs[eth], SOFTPLUS_V_EPS) * cX[i,j] / cs[IDX_E]

        ode_O2[i=1:nfe, j=1:ncp],
            cdot[IDX_O2,i,j] == -MW_O2 * (-v[o2,i] * vs[o2]) * cX[i,j] / cs[IDX_O2]

        ode_Prot[i=1:nfe, j=1:ncp],
            cdot[IDX_PROT,i,j] ==
                (softplus_jl(v[obj,i] * vs[obj], SOFTPLUS_V_EPS) * cX[i,j] * PROT_CONTENT_0
                 - cProt[i,j] * K_DEATH
                 + XA_FRACTION * cX[i,j] * softplus_jl(v[IDX_PROT_RXN,i] * vs[IDX_PROT_RXN], SOFTPLUS_V_EPS)
                 - TURNOVER_LAMBDA * cProt[i,j]) / cs[IDX_PROT]

        ode_Carb[i=1:nfe, j=1:ncp],
            cdot[IDX_CARB,i,j] ==
                (softplus_jl(v[obj,i] * vs[obj], SOFTPLUS_V_EPS) * cX[i,j] * CARB_CONTENT_0
                 - cCarb[i,j] * K_DEATH) / cs[IDX_CARB]

        # dynamic sugar / kinetic-N uptake caps at FE endpoint
        dyn_upt[k=1:N_UP, i=1:nfe],
            -v[UPTAKE_IDXS[k],i] * vs[UPTAKE_IDXS[k]] - phase_g_fe[i] * L_uptake[k,i,ncp] <= 0

        dyn_prod[k=1:N_PROD, i=1:nfe],
            v[PRODUCT_IDXS[k],i] * vs[PRODUCT_IDXS[k]] - L_product_cap[k,i] <= 0

        atpm_floor[i=1:nfe],
            -v[IDX_ATPM,i] * vs[IDX_ATPM] + ATPM_floor_rhs[i] <= 0

        # complementarity base bounds
        FO1[mc=1:nv, i=1:nfe],
            FO_L[mc,i] == (v[mc,i] * vs[mc] - vlb[mc]) * alpha_L[mc,i]

        FO2[mc=1:nv, i=1:nfe],
            FO_U[mc,i] == (v[mc,i] * vs[mc] - vub[mc]) * alpha_U[mc,i]

        FO3[k=1:N_UP, i=1:nfe],
            FO_upt[k,i] == (-v[UPTAKE_IDXS[k],i] * vs[UPTAKE_IDXS[k]] - phase_g_fe[i] * L_uptake[k,i,ncp]) * alpha_upt[k,i]

        FO4[k=1:N_PROD, i=1:nfe],
            FO_prod[k,i] == (v[PRODUCT_IDXS[k],i] * vs[PRODUCT_IDXS[k]] - L_product_cap[k,i]) * alpha_prod[k,i]

        FO_atpm_c[i=1:nfe],
            FO_atpm[i] == (-v[IDX_ATPM,i] * vs[IDX_ATPM] + ATPM_floor_rhs[i]) * alpha_atpm_floor[i]
    end)

    if N_AA > 0
        @NLconstraints(m, begin
            aa_uptake_cap[a=1:N_AA, i=1:nfe],
                -v[AA_UPTAKE_IDXS[a],i] * vs[AA_UPTAKE_IDXS[a]] - Q_AA_cap[a,i] <= 0

            FO_aa_c[a=1:N_AA, i=1:nfe],
                FO_aa[a,i] == (-v[AA_UPTAKE_IDXS[a],i] * vs[AA_UPTAKE_IDXS[a]] - Q_AA_cap[a,i]) * alpha_aa[a,i]

            aa_state_dyn[a=1:N_AA, i=1:nfe, j=1:ncp],
                cdot[AA_STATE_IDXS[a],i,j] == phase_g[i,j] * (v[AA_UPTAKE_IDXS[a],i] * vs[AA_UPTAKE_IDXS[a]]) * cX[i,j] / cs[AA_STATE_IDXS[a]]
        end)
    end

    if N_ARO > 0
        @NLconstraint(m, aroma_state_dyn[a=1:N_ARO, i=1:nfe, j=1:ncp],
            cdot[AROMA_STATE_IDXS[a],i,j] ==
                AROMA_MW_VEC[a] * softplus_jl(v[AROMA_RXN_IDXS[a],i] * vs[AROMA_RXN_IDXS[a]], SOFTPLUS_V_EPS) * cX[i,j] / cs[AROMA_STATE_IDXS[a]])
    end

    if N_PAIR > 0
        @NLconstraints(m, begin
            pair_c[p=1:N_PAIR, i=1:nfe],
                PAIRWISE_PHI[p] * v[PAIRWISE_ALCOHOL_IDXS[p],i] * vs[PAIRWISE_ALCOHOL_IDXS[p]] -
                v[PAIRWISE_ESTER_IDXS[p],i] * vs[PAIRWISE_ESTER_IDXS[p]] <= 0

            FO_pair_c[p=1:N_PAIR, i=1:nfe],
                FO_pair[p,i] ==
                    (PAIRWISE_PHI[p] * v[PAIRWISE_ALCOHOL_IDXS[p],i] * vs[PAIRWISE_ALCOHOL_IDXS[p]] -
                     v[PAIRWISE_ESTER_IDXS[p],i] * vs[PAIRWISE_ESTER_IDXS[p]]) * alpha_pair[p,i]
        end)
    end

    if HAS_EA_SOFT
        @NLconstraints(m, begin
            ea_c[i=1:nfe],
                PHI_ETHYL_ACETATE_STATIC * v[EA_SOFT_ALCOHOL_IDX,i] * vs[EA_SOFT_ALCOHOL_IDX] -
                v[EA_SOFT_ESTER_IDX,i] * vs[EA_SOFT_ESTER_IDX] <= 0

            FO_ea_c[i=1:nfe],
                FO_ea[i] ==
                    (PHI_ETHYL_ACETATE_STATIC * v[EA_SOFT_ALCOHOL_IDX,i] * vs[EA_SOFT_ALCOHOL_IDX] -
                     v[EA_SOFT_ESTER_IDX,i] * vs[EA_SOFT_ESTER_IDX]) * alpha_ea[i]
        end)
    end

    if !HAS_EA_SOFT
        @constraints(m, begin
            ea_alpha_off[i=1:nfe], alpha_ea[i] == 0.0
            ea_fo_off[i=1:nfe], FO_ea[i] == 0.0
        end)
    end

    optimize!(m)

    term = termination_status(m)
    prim = primal_status(m)
    dual = dual_status(m)

    println("\nStatus = ", term)
    println("Primal = ", prim, " | Dual = ", dual)

    cStar_raw = value.(c)
    cStar = similar(cStar_raw)
    for s in 1:nc, i in 1:nfe, j in 1:ncp
        cStar[s,i,j] = cStar_raw[s,i,j] * cs[s]
    end
    vStar_raw = value.(v)
    vStar = similar(vStar_raw)
    for rx in 1:nv, i in 1:nfe
        vStar[rx,i] = vStar_raw[rx,i] * vs[rx]
    end
    cdStar  = value.(cdot)
    lStar   = value.(lambda_)
    alLStar = value.(alpha_L)
    alUStar = value.(alpha_U)
    auStar  = value.(alpha_upt)
    apStar  = value.(alpha_prod)
    aaaStar = value.(alpha_aa)
    aprStar = value.(alpha_pair)
    aeaStar = value.(alpha_ea)
    aatpmStar = value.(alpha_atpm_floor)
    hStar   = value.(hv)

    diagStar = (
        solver = (term = term, primal = prim, dual = dual),
        phase_growth = value.(phase_g_fe),
        gam_extra = GAM_EXTRA_FE,
        atpm_rhs = value.(ATPM_floor_rhs),
        L_upt = N_UP > 0 ? value.(L_uptake[:, :, ncp]) : zeros(0, nfe),
        L_prod = N_PROD > 0 ? value.(L_product[:, :, ncp]) : zeros(0, nfe),
        L_prod_cap = N_PROD > 0 ? value.(L_product_cap) : zeros(0, nfe),
        Q_AA_cap = N_AA > 0 ? value.(Q_AA_cap) : zeros(0, nfe),
        uptake_violation = N_UP > 0 ? maximum([-vStar[UPTAKE_IDXS[k],i] - value(L_uptake[k,i,ncp]) for k in 1:N_UP, i in 1:nfe]) : 0.0,
        product_violation = N_PROD > 0 ? maximum([vStar[PRODUCT_IDXS[k],i] - value(L_product_cap[k,i]) for k in 1:N_PROD, i in 1:nfe]) : 0.0,
        pairwise_active = N_PAIR > 0,
        ethyl_acetate_soft = HAS_EA_SOFT,
    )

    return cStar, vStar, cdStar, lStar, alLStar, alUStar, auStar, apStar, aaaStar, aprStar, aeaStar, aatpmStar, hStar, diagStar
end
