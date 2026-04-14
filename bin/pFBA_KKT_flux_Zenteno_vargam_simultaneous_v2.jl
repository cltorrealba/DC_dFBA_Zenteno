# ======================================================================================
# pFBA_KKT_flux_Zenteno_vargam_simultaneous_v2.jl
# --------------------------------------------------------------------------------------
# Version v2: homologa exactamente la formulacion secuencial Radau IIA de la celda 46
# del notebook 1 (variante balanceada _run_dfba_colloc_balanced). Cambios clave respecto
# al archivo base:
#
#   1) Nuevo estado IDX_NREC (pool de nitrogeno reciclado por turnover de proteina).
#      - En fase BIOMASS, ecuacion de N_rec:  d(N_rec)/dt = Y_N_FROM_PROT * lambda * cProt
#      - En fase TURNOVER_ATPM, ademas se drena N_rec para fuelear el uptake de AA:
#           d(N_rec)/dt -= (sum_a v_aa_a) * MW_N * cX
#        El uptake de AA en turnover esta limitado por
#           q_AA_a <= (K_NREC_UPTAKE * N_rec / cX) * AA_ALPHA[a] / MW_N
#   2) Conmutacion de fase sigue siendo suave (sigmoidal) en cN_total = cN_free + cNaa,
#      igual que el archivo base, pero el cap de uptake de AA en turnover usa N_rec.
#   3) Fase BIOMASS: ATPM queda liberado (BAL_ATPM_GROWTH_LB, BAL_ATPM_GROWTH_UB) y se
#      mantiene el bloqueo mu >= BAL_GROWTH_LOCK_FRAC * mu* implicito via el peso de
#      crecimiento en la funcion de objetivo dinamica (d_phase) y el piso atpm.
#   4) GAM variable externa via GAM_EXTRA_FE (exportada desde NB1, identica al base).
#   5) Se agrega el kwarg opcional `seed::Union{Nothing,Dict}=nothing` para sembrar
#      warm-start en CADA variable de decision del MPCC (c, cdot, v, lambda_,
#      alpha_*, FO_*, hv). Si una clave del dict falta se omite con @info.
#   6) Se conserva la estructura KKT / MPCC completa (complementariedades FO_*,
#      stationarity, product/uptake/AA/pair/EA/atpm caps) del archivo base.
#
# Globales esperadas (ademas de las del archivo base):
#   IDX_NREC, Y_N_FROM_PROT, K_NREC_UPTAKE, BAL_ATPM_GROWTH_LB, BAL_ATPM_GROWTH_UB
# ======================================================================================

using JuMP
using Ipopt
using LinearAlgebra
using Statistics
using ForwardDiff
import MathOptInterface as MOI

_softplus(x, ϵ) = 0.5 * (x + sqrt(x * x + ϵ * ϵ))
_sigmoid(x) = 1.0 / (1.0 + exp(-x))

# ──────────────────────────────────────────────────────────────────────────────
# Derivadas analíticas / ForwardDiff para JuMP.register  (habilita eval_h)
# ──────────────────────────────────────────────────────────────────────────────

# ── sigmoid: f(x), f'(x), f''(x) ──
_sigmoid_d1(x) = (s = _sigmoid(x); s * (1.0 - s))
_sigmoid_d2(x) = (s = _sigmoid(x); ds = s * (1.0 - s); ds * (1.0 - 2.0 * s))

# ── softplus(x, ε): grad 2-vector, Hess 2×2 (lower-triangle) ──
function _softplus_grad(g, x, ε)
    r = sqrt(x * x + ε * ε)
    g[1] = 0.5 * (1.0 + x / r)
    g[2] = 0.5 * ε / r
end
function _softplus_hess(H, x, ε)
    r3 = (x * x + ε * ε)^1.5
    H[1,1] =  0.5 * ε * ε / r3
    H[2,1] = -0.5 * x * ε / r3
    H[2,2] =  0.5 * x * x / r3
end

# ── dynamic_temperature: f(t), f'(t), f''(t) ──
#    Definidas como closures que capturan T_BASE, T_STEPS, T_DELTAS, T_STEEP desde Main
function _dtemp_d1(t)
    val = 0.0
    for i in eachindex(Main.T_STEPS)
        σ = 1.0 / (1.0 + exp(-Main.T_STEEP * (t - Main.T_STEPS[i])))
        val += Main.T_DELTAS[i] * Main.T_STEEP * σ * (1.0 - σ)
    end
    return val
end
function _dtemp_d2(t)
    val = 0.0
    for i in eachindex(Main.T_STEPS)
        σ = 1.0 / (1.0 + exp(-Main.T_STEEP * (t - Main.T_STEPS[i])))
        ds = σ * (1.0 - σ)
        val += Main.T_DELTAS[i] * Main.T_STEEP^2 * ds * (1.0 - 2.0 * σ)
    end
    return val
end

# ── death_rate_T(E, T): gradiente y hessiana via ForwardDiff ──
function _death_rate_T_vec(x::AbstractVector)
    E     = x[1]
    T_val = x[2]
    _R    = Main.R
    _Kd0  = Main.Kd0_nom
    Td   = -0.0001 * E^3 + 0.0049 * E^2 - 0.1279 * E + 315.89
    s    = 0.5 * (1.0 + tanh(0.5 * (T_val - Td)))
    base = _Kd0 * exp(0.0415 * E + (130000.0 * (T_val - 305.65)) / (305.65 * _R * T_val))
    return base * s
end
function _death_grad(g, E, T_val)
    gv = ForwardDiff.gradient(_death_rate_T_vec, [E, T_val])
    g[1] = gv[1]; g[2] = gv[2]
end
function _death_hess(H, E, T_val)
    Hv = ForwardDiff.hessian(_death_rate_T_vec, [E, T_val])
    H[1,1] = Hv[1,1]; H[2,1] = Hv[2,1]; H[2,2] = Hv[2,2]
end

# ── smooth_injection(t, t_shot, dose, width): grad 4-vec, Hess 4×4 via ForwardDiff ──
function _injection_vec(x::AbstractVector)
    t, t_shot, dose, width = x[1], x[2], x[3], x[4]
    _SQRT_2PI = Main.SQRT_2PI
    return (dose / (width * _SQRT_2PI)) * exp(-0.5 * ((t - t_shot) / width)^2)
end
function _injection_grad(g, t, t_shot, dose, width)
    gv = ForwardDiff.gradient(_injection_vec, [t, t_shot, dose, width])
    for i in 1:4; g[i] = gv[i]; end
end
function _injection_hess(H, t, t_shot, dose, width)
    Hv = ForwardDiff.hessian(_injection_vec, [t, t_shot, dose, width])
    for i in 1:4, j in 1:i
        H[i,j] = Hv[i,j]
    end
end

# --------------------------------------------------------------------------------------
# Helper: setear start values de forma segura desde un dict
# --------------------------------------------------------------------------------------
function _seed_apply!(var, arr, name::Symbol)
    try
        sz_v = size(var)
        sz_a = size(arr)
        if sz_v != sz_a
            @info "[seed] shape mismatch en $(name): var=$(sz_v) seed=$(sz_a); se omite."
            return
        end
        for I in eachindex(var)
            val = arr[I]
            if isfinite(val)
                set_start_value(var[I], val)
            end
        end
    catch err
        @info "[seed] no se pudo aplicar $(name): $(err)"
    end
end

function _seed_get(seed::Dict, key::Symbol)
    if haskey(seed, key)
        return seed[key]
    elseif haskey(seed, String(key))
        return seed[String(key)]
    else
        return nothing
    end
end

# --------------------------------------------------------------------------------------
# Funcion principal
# --------------------------------------------------------------------------------------
function pFBA_KKT_flux_Zenteno_vargam_simultaneous_v2(
    c0;
    eps_flux::Float64 = 0.0,
    apply_product_caps::Bool = false,
    seed::Union{Nothing,Dict} = nothing,
    eps_schedule::Vector{Float64} = [1.0, 1e-1, 1e-2, 1e-3, 1e-4],
    verbose_homotopy::Bool = true,
)
    @assert length(c0) == nc "c0 debe tener longitud nc."
    @assert nc >= 9 "Este modelo v2 asume al menos 9 estados base (incluye IDX_NREC)."

    # ---------- metadata sizes ----------
    N_UP   = length(UPTAKE_IDXS)
    N_PROD = length(PRODUCT_IDXS)
    N_AA   = length(AA_UPTAKE_IDXS)
    N_ARO  = length(AROMA_RXN_IDXS)
    N_PAIR = length(PAIRWISE_ESTER_IDXS)
    HAS_EA_SOFT = (EA_SOFT_ESTER_IDX > 0) && (EA_SOFT_ALCOHOL_IDX > 0) && (PHI_ETHYL_ACETATE_STATIC > 0.0)

    # ---------- sparse structural maps (critical for IPOPT memory/speed) ----------
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

    # ---------- Constantes opcionales (defaults si no estan definidas globalmente) ----------
    _Y_N_FROM_PROT      = isdefined(Main, :Y_N_FROM_PROT)      ? Y_N_FROM_PROT      : 0.16
    _K_NREC_UPTAKE      = isdefined(Main, :K_NREC_UPTAKE)      ? K_NREC_UPTAKE      : 0.05
    _BAL_ATPM_GROW_LB   = isdefined(Main, :BAL_ATPM_GROWTH_LB) ? BAL_ATPM_GROWTH_LB : 0.7
    _BAL_ATPM_GROW_UB   = isdefined(Main, :BAL_ATPM_GROWTH_UB) ? BAL_ATPM_GROWTH_UB : 1000.0

    # ---------- Radau IIA ----------
    colmat = [0.19681547722366  -0.06553542585020   0.02377097434822;
              0.39442431473909   0.29207341166523  -0.04154875212600;
              0.37640306270047   0.51248582618842   0.11111111111111]
    radau  = [0.15505102572168, 0.64494897427832, 1.0]

    # ---------- Model ----------
    m = Model(optimizer_with_attributes(
        Ipopt.Optimizer,
        "linear_solver" => IPOPT_LINEAR_SOLVER,
        "hsllib"        => get(ENV, "IPOPT_HSLLIB", ""),
        "print_level"   => IPOPT_PRINT_LEVEL,
        "tol"           => IPOPT_TOL,
        "acceptable_tol"=> IPOPT_ACCEPTABLE_TOL,
        "acceptable_iter" => IPOPT_ACCEPTABLE_ITER,
        "max_iter"      => IPOPT_MAX_ITER,
        # mu_strategy: monotone con mu_init alto es mas robusto frente a inf_pr inicial
        # alto que adaptive (que se asusta y entra a restoration de inmediato).
        "mu_strategy"   => "monotone",
        "mu_init"       => 1e-1,
        "mu_min"        => 1e-9,
        "constr_viol_tol" => IPOPT_CONSTR_VIOL_TOL,
        "compl_inf_tol" => IPOPT_COMPL_INF_TOL,
        "bound_relax_factor" => 1e-8,
        "hessian_approximation" => "exact",
        "nlp_scaling_method"    => "gradient-based",
        "alpha_for_y"           => "min",
        "recalc_y"              => "yes",
        # Tolerar infactibilidad inicial alta sin saltar a restoration:
        "expect_infeasible_problem" => "yes",
        "expect_infeasible_problem_ctol" => 1e-2,
        "expect_infeasible_problem_ytol" => 1e8,
        "required_infeasibility_reduction" => 0.5,
        "soft_resto_pderror_reduction_factor" => 0.999,
        "warm_start_init_point" => "yes",
        "warm_start_bound_push" => 1e-6,
        "warm_start_mult_bound_push" => 1e-6,
        "warm_start_bound_frac" => 1e-6,
        "warm_start_slack_bound_push" => 1e-6,
        "warm_start_slack_bound_frac" => 1e-6,
        "print_user_options"    => "yes",
    ))

    # output_file se agrega solo si la ruta es ASCII valida (IPOPT rechaza non-ASCII)
    if isdefined(Main, :IPOPT_OUTPUT_FILE) && !isempty(IPOPT_OUTPUT_FILE) && isascii(IPOPT_OUTPUT_FILE)
        set_optimizer_attribute(m, "output_file", IPOPT_OUTPUT_FILE)
        @info "[IPOPT] output_file = $IPOPT_OUTPUT_FILE"
    end

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

    # ---------- starts default (sobrescritos por seed si viene) ----------
    for i in 1:nfe
        set_start_value(hv[i], hm[i])
        for j in 1:ncp, s in 1:nc
            set_start_value(c[s, i, j], c0[s] / cs[s])
        end
        for rx in 1:nv
            set_start_value(v[rx, i], 0.0)
        end
    end

    # ---------- aplicar seed completo si fue provisto ----------
    if seed !== nothing
        @info "[seed] Aplicando warm-start desde dict a todas las variables del MPCC..."
        _map = Dict(
            :c=>c, :cdot=>cdot, :v=>v, :lambda_=>lambda_,
            :alpha_L=>alpha_L, :alpha_U=>alpha_U,
            :alpha_upt=>alpha_upt, :alpha_prod=>alpha_prod,
            :alpha_aa=>alpha_aa, :alpha_pair=>alpha_pair,
            :alpha_ea=>alpha_ea, :alpha_atpm_floor=>alpha_atpm_floor,
            :FO_L=>FO_L, :FO_U=>FO_U,
            :FO_upt=>FO_upt, :FO_prod=>FO_prod,
            :FO_aa=>FO_aa, :FO_pair=>FO_pair,
            :FO_ea=>FO_ea, :FO_atpm=>FO_atpm,
            :hv=>hv,
        )
        for (k, jv) in _map
            arr = _seed_get(seed, k)
            if arr === nothing
                @info "[seed] clave ausente: $(k) (se omite)"
            else
                _seed_apply!(jv, arr, k)
            end
        end
    end

    c0s = [c0[i] / cs[i] for i in 1:nc]

    # ---------- register helper functions (con derivadas explicitas → Hessian exacto) ----------
    JuMP.register(m, :dynamic_temperature, 1, dynamic_temperature, _dtemp_d1, _dtemp_d2)
    JuMP.register(m, :death_rate_T, 2, death_rate_T, _death_grad, _death_hess)
    JuMP.register(m, :smooth_injection, 4, smooth_injection, _injection_grad, _injection_hess)
    JuMP.register(m, :softplus_jl, 2, _softplus, _softplus_grad, _softplus_hess)
    JuMP.register(m, :sigmoid_jl, 1, _sigmoid, _sigmoid_d1, _sigmoid_d2)

    # ---------- mapa temporal ----------
    @NLexpression(m, t_fe[i=1:nfe], sum(hv[k] for k in 1:i-1))
    @NLexpression(m, t_loc[i=1:nfe, j=1:ncp], t_fe[i] + radau[j] * hv[i])

    # ---------- alias escalado->fisico ----------
    @NLexpressions(m, begin
        cX[i=1:nfe, j=1:ncp], c[IDX_X, i, j] * cs[IDX_X]
        cNf[i=1:nfe, j=1:ncp], c[IDX_NFREE, i, j] * cs[IDX_NFREE]
        cG[i=1:nfe, j=1:ncp], c[IDX_G, i, j] * cs[IDX_G]
        cF[i=1:nfe, j=1:ncp], c[IDX_F, i, j] * cs[IDX_F]
        cE[i=1:nfe, j=1:ncp], c[IDX_E, i, j] * cs[IDX_E]
        cO2[i=1:nfe, j=1:ncp], c[IDX_O2, i, j] * cs[IDX_O2]
        cProt[i=1:nfe, j=1:ncp], c[IDX_PROT, i, j] * cs[IDX_PROT]
        cCarb[i=1:nfe, j=1:ncp], c[IDX_CARB, i, j] * cs[IDX_CARB]
        cNrec[i=1:nfe, j=1:ncp], c[IDX_NREC, i, j] * cs[IDX_NREC]
        cNaa[i=1:nfe, j=1:ncp], MW_N * sum(c[AA_STATE_IDXS[a], i, j] * cs[AA_STATE_IDXS[a]] for a in 1:N_AA)
        cNtot[i=1:nfe, j=1:ncp], cNf[i,j] + cNaa[i,j]
        phase_g[i=1:nfe, j=1:ncp], sigmoid_jl((cNtot[i,j] - N_TOTAL_DEPLETION_THRESHOLD) / PHASE_SMOOTH_EPS)
        phase_t[i=1:nfe, j=1:ncp], 1.0 - phase_g[i,j]
    end)

    @NLexpression(m, phase_g_fe[i=1:nfe], phase_g[i, ncp])
    @NLexpression(m, phase_t_fe[i=1:nfe], 1.0 - phase_g_fe[i])

    # ---------- envolventes cineticas / temperatura (igual al base) ----------
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

    # ---------- limites cineticos agregados ----------
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

    # ---------- Caps de AA: fase crecimiento vs turnover (NUEVO: usa N_rec) ----------
    # Crecimiento:  q_AA <= K_AA_UPTAKE_GROWTH * cAA / cX
    # Turnover  :  q_AA <= (K_NREC_UPTAKE * N_rec / cX) * AA_ALPHA[a] / MW_N
    if N_AA > 0
        @NLexpression(m, Q_AA_growth[a=1:N_AA, i=1:nfe],
            K_AA_UPTAKE_GROWTH * (c[AA_STATE_IDXS[a], i, ncp] * cs[AA_STATE_IDXS[a]]) / (cX[i,ncp] + EPS))
        @NLexpression(m, Q_AA_turn[a=1:N_AA, i=1:nfe],
            (_K_NREC_UPTAKE * cNrec[i,ncp] / (cX[i,ncp] + EPS)) * AA_ALPHA_VEC[a] / MW_N)
        @NLexpression(m, Q_AA_cap[a=1:N_AA, i=1:nfe],
            phase_g_fe[i] * Q_AA_growth[a,i] + phase_t_fe[i] * Q_AA_turn[a,i])
    end

    # ---------- ATPM floor con GAM variable ----------
    # NOTA: en la variante balanceada la fase BIOMASS libera ATPM; se modela como piso
    # _BAL_ATPM_GROW_LB durante crecimiento (no como igualdad).
    @NLexpression(m, ATPM_floor_rhs[i=1:nfe],
        phase_t_fe[i] * ATPM_LB_NO_GROWTH +
        phase_g_fe[i] * _BAL_ATPM_GROW_LB +
        GAM_EXTRA_FE[i] * softplus_jl(v[obj,i] * vs[obj], SOFTPLUS_V_EPS))

    # ---------- Coeficientes dinamicos de objetivo ----------
    @NLexpression(m, d_phase[mc=1:nv, i=1:nfe],
        phase_g_fe[i] * D_GROWTH[mc] + phase_t_fe[i] * D_TURNOVER[mc] + phase_g_fe[i] * D_AAUP[mc])

    # ---------- Objetivo externo (Scholtes MPCC: solo regularizacion pFBA) ----------
    # Sin penalty sobre FO_*; complementariedades se manejan via bounds |FO_*| <= eps
    # relajadas en homotopia (Scholtes).
    @NLobjective(m, Min,
        sum(
            Q_REG * sum((v[mc,i] * vs[mc])^2 for mc in 1:nv) +
            FLUX_SMOOTH_WEIGHT * sum((v[mc,i] * vs[mc])^2 for mc in 1:nv)
        for i in 1:nfe)
    )

    # ---------- Restricciones lineales ----------
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

    # ---------- Sc (sparse) ----------
    Sc = Matrix{ConstraintRef}(undef, nm, nfe)
    for mc in 1:nm
        cols = S_ROW_COLS[mc]
        vals = S_ROW_VALS[mc]
        for i in 1:nfe
            Sc[mc,i] = @constraint(m, sum(vals[t] * v[cols[t], i] * vs[cols[t]] for t in eachindex(cols)) == 0)
        end
    end

    # ---------- Stationarity Lagrangiana (sparse) ----------
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

    # ---------- ODEs dinamicas (mirror de cell 46 balanced) ----------
    @NLconstraints(m, begin
        # Biomasa: crece solo en fase BIOMASS (gated por phase_g)
        ode_X[i=1:nfe, j=1:ncp],
            cdot[IDX_X,i,j] == ((phase_g[i,j] * softplus_jl(v[obj,i] * vs[obj], SOFTPLUS_V_EPS) - Kd_j[i,j]) * cX[i,j]) / cs[IDX_X]

        # N libre: solo consumido en fase BIOMASS; incluye injection (DAP dosing)
        ode_Nfree[i=1:nfe, j=1:ncp],
            cdot[IDX_NFREE,i,j] ==
                (- MW_N * phase_g[i,j] *
                    sum(IS_NIT[k] * N_atoms_vec[UPTAKE_IDXS[k]] *
                        (-v[UPTAKE_IDXS[k],i] * vs[UPTAKE_IDXS[k]]) for k in 1:N_UP) * cX[i,j]
                 + injection_rate[i,j]) / cs[IDX_NFREE]

        # Glucosa, Fructosa: consumo en ambas fases
        ode_G[i=1:nfe, j=1:ncp],
            cdot[IDX_G,i,j] == -MW_GLU * (-v[glu,i] * vs[glu]) * cX[i,j] / cs[IDX_G]

        ode_F[i=1:nfe, j=1:ncp],
            cdot[IDX_F,i,j] == -MW_FRU * (-v[fru,i] * vs[fru]) * cX[i,j] / cs[IDX_F]

        ode_E[i=1:nfe, j=1:ncp],
            cdot[IDX_E,i,j] == MW_ETH * softplus_jl(v[eth,i] * vs[eth], SOFTPLUS_V_EPS) * cX[i,j] / cs[IDX_E]

        ode_O2[i=1:nfe, j=1:ncp],
            cdot[IDX_O2,i,j] == -MW_O2 * (-v[o2,i] * vs[o2]) * cX[i,j] / cs[IDX_O2]

        # Proteina con turnover (Henriques 2023): sintesis + Xa*v_prot - death - turnover
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

        # ================= NUEVO ESTADO: N_rec =================
        # Entrada: turnover de proteina libera N reciclado (Y_N_FROM_PROT * lambda * cProt)
        # Salida : solo en fase turnover, fuelear uptake de AA:
        #           sum_a (-v_aa_a) * MW_N * cX
        ode_Nrec[i=1:nfe, j=1:ncp],
            cdot[IDX_NREC,i,j] ==
                (_Y_N_FROM_PROT * TURNOVER_LAMBDA * cProt[i,j]
                 - phase_t[i,j] * MW_N * cX[i,j] *
                    sum((-v[AA_UPTAKE_IDXS[a],i] * vs[AA_UPTAKE_IDXS[a]]) for a in 1:N_AA)
                ) / cs[IDX_NREC]

        # ---------- Caps dinamicos (uptake, producto, ATPM) ----------
        dyn_upt[k=1:N_UP, i=1:nfe],
            -v[UPTAKE_IDXS[k],i] * vs[UPTAKE_IDXS[k]] - phase_g_fe[i] * L_uptake[k,i,ncp] <= 0

        dyn_prod[k=1:N_PROD, i=1:nfe],
            v[PRODUCT_IDXS[k],i] * vs[PRODUCT_IDXS[k]] - L_product_cap[k,i] <= 0

        atpm_floor[i=1:nfe],
            -v[IDX_ATPM,i] * vs[IDX_ATPM] + ATPM_floor_rhs[i] <= 0

        # ---------- Complementariedades base ----------
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

    # ---------- AA dynamics ----------
    if N_AA > 0
        @NLconstraints(m, begin
            aa_uptake_cap[a=1:N_AA, i=1:nfe],
                -v[AA_UPTAKE_IDXS[a],i] * vs[AA_UPTAKE_IDXS[a]] - Q_AA_cap[a,i] <= 0

            FO_aa_c[a=1:N_AA, i=1:nfe],
                FO_aa[a,i] == (-v[AA_UPTAKE_IDXS[a],i] * vs[AA_UPTAKE_IDXS[a]] - Q_AA_cap[a,i]) * alpha_aa[a,i]

            # En cell 46: AA state varia en ambas fases (crecimiento consume,
            # turnover tambien consume via N_rec). Dejamos el acoplamiento completo:
            aa_state_dyn[a=1:N_AA, i=1:nfe, j=1:ncp],
                cdot[AA_STATE_IDXS[a],i,j] == (v[AA_UPTAKE_IDXS[a],i] * vs[AA_UPTAKE_IDXS[a]]) * cX[i,j] / cs[AA_STATE_IDXS[a]]
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

    # ---------- Scholtes relaxation: |FO_*| <= eps ----------
    # Los refs se actualizan en el loop de homotopia via set_normalized_rhs.
    eps_init = eps_schedule[1]
    @constraints(m, begin
        sch_FO_L_u[mc=1:nv, i=1:nfe],  FO_L[mc,i] <=  eps_init
        sch_FO_L_l[mc=1:nv, i=1:nfe],  FO_L[mc,i] >= -eps_init
        sch_FO_U_u[mc=1:nv, i=1:nfe],  FO_U[mc,i] <=  eps_init
        sch_FO_U_l[mc=1:nv, i=1:nfe],  FO_U[mc,i] >= -eps_init
        sch_FO_atpm_u[i=1:nfe], FO_atpm[i] <=  eps_init
        sch_FO_atpm_l[i=1:nfe], FO_atpm[i] >= -eps_init
        sch_FO_ea_u[i=1:nfe],   FO_ea[i]   <=  eps_init
        sch_FO_ea_l[i=1:nfe],   FO_ea[i]   >= -eps_init
    end)
    if N_UP > 0
        @constraints(m, begin
            sch_FO_upt_u[k=1:N_UP, i=1:nfe], FO_upt[k,i] <=  eps_init
            sch_FO_upt_l[k=1:N_UP, i=1:nfe], FO_upt[k,i] >= -eps_init
        end)
    end
    if N_PROD > 0
        @constraints(m, begin
            sch_FO_prod_u[k=1:N_PROD, i=1:nfe], FO_prod[k,i] <=  eps_init
            sch_FO_prod_l[k=1:N_PROD, i=1:nfe], FO_prod[k,i] >= -eps_init
        end)
    end
    if N_AA > 0
        @constraints(m, begin
            sch_FO_aa_u[a=1:N_AA, i=1:nfe], FO_aa[a,i] <=  eps_init
            sch_FO_aa_l[a=1:N_AA, i=1:nfe], FO_aa[a,i] >= -eps_init
        end)
    end
    if N_PAIR > 0
        @constraints(m, begin
            sch_FO_pair_u[p=1:N_PAIR, i=1:nfe], FO_pair[p,i] <=  eps_init
            sch_FO_pair_l[p=1:N_PAIR, i=1:nfe], FO_pair[p,i] >= -eps_init
        end)
    end

    function _update_eps!(eps::Float64)
        for mc in 1:nv, i in 1:nfe
            set_normalized_rhs(sch_FO_L_u[mc,i],  eps)
            set_normalized_rhs(sch_FO_L_l[mc,i], -eps)
            set_normalized_rhs(sch_FO_U_u[mc,i],  eps)
            set_normalized_rhs(sch_FO_U_l[mc,i], -eps)
        end
        for i in 1:nfe
            set_normalized_rhs(sch_FO_atpm_u[i],  eps)
            set_normalized_rhs(sch_FO_atpm_l[i], -eps)
            set_normalized_rhs(sch_FO_ea_u[i],    eps)
            set_normalized_rhs(sch_FO_ea_l[i],   -eps)
        end
        if N_UP > 0
            for k in 1:N_UP, i in 1:nfe
                set_normalized_rhs(sch_FO_upt_u[k,i],  eps)
                set_normalized_rhs(sch_FO_upt_l[k,i], -eps)
            end
        end
        if N_PROD > 0
            for k in 1:N_PROD, i in 1:nfe
                set_normalized_rhs(sch_FO_prod_u[k,i],  eps)
                set_normalized_rhs(sch_FO_prod_l[k,i], -eps)
            end
        end
        if N_AA > 0
            for a in 1:N_AA, i in 1:nfe
                set_normalized_rhs(sch_FO_aa_u[a,i],  eps)
                set_normalized_rhs(sch_FO_aa_l[a,i], -eps)
            end
        end
        if N_PAIR > 0
            for p in 1:N_PAIR, i in 1:nfe
                set_normalized_rhs(sch_FO_pair_u[p,i],  eps)
                set_normalized_rhs(sch_FO_pair_l[p,i], -eps)
            end
        end
        return nothing
    end

    # ---------- Homotopia Scholtes ----------
    homotopy_log = NamedTuple[]
    term = MOI.OPTIMIZE_NOT_CALLED
    prim = MOI.NO_SOLUTION
    dual = MOI.NO_SOLUTION
    for (k_eps, eps_k) in enumerate(eps_schedule)
        if k_eps > 1
            _update_eps!(eps_k)
        end
        if verbose_homotopy
            println("\n================================================================")
            println("[homotopy] iter $k_eps / $(length(eps_schedule))  eps = $eps_k")
            println("================================================================")
        end
        optimize!(m)
        term = termination_status(m)
        prim = primal_status(m)
        dual = dual_status(m)
        obj_k = try objective_value(m) catch; NaN end
        push!(homotopy_log, (iter=k_eps, eps=eps_k, term=term, prim=prim, dual=dual, obj=obj_k))
        if verbose_homotopy
            println("[homotopy] iter $k_eps eps=$eps_k status=$term primal=$prim obj=$obj_k")
        end
        if term in (MOI.INFEASIBLE, MOI.LOCALLY_INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED)
            @warn "[homotopy] infactibilidad detectada en eps=$eps_k. Deteniendo homotopia."
            break
        end
        if prim != MOI.FEASIBLE_POINT && prim != MOI.NEARLY_FEASIBLE_POINT
            @warn "[homotopy] primal no factible en eps=$eps_k ($prim). Deteniendo para inspeccion."
            break
        end
    end

    println("\n[homotopy] Resumen:")
    for row in homotopy_log
        println("  eps=$(row.eps)  term=$(row.term)  primal=$(row.prim)  obj=$(row.obj)")
    end
    println("\nStatus final = ", term)
    println("Primal = ", prim, " | Dual = ", dual)

    # ---------- Desescalado ----------
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
    cdStar    = value.(cdot)
    lStar     = value.(lambda_)
    alLStar   = value.(alpha_L)
    alUStar   = value.(alpha_U)
    auStar    = value.(alpha_upt)
    apStar    = value.(alpha_prod)
    aaaStar   = value.(alpha_aa)
    aprStar   = value.(alpha_pair)
    aeaStar   = value.(alpha_ea)
    aatpmStar = value.(alpha_atpm_floor)
    hStar     = value.(hv)

    diagStar = (
        solver = (term = term, primal = prim, dual = dual),
        phase_growth = value.(phase_g_fe),
        gam_extra = GAM_EXTRA_FE,
        atpm_rhs = value.(ATPM_floor_rhs),
        L_upt = N_UP > 0 ? value.(L_uptake[:, :, ncp]) : zeros(0, nfe),
        L_prod = N_PROD > 0 ? value.(L_product[:, :, ncp]) : zeros(0, nfe),
        L_prod_cap = N_PROD > 0 ? value.(L_product_cap) : zeros(0, nfe),
        Q_AA_cap = N_AA > 0 ? value.(Q_AA_cap) : zeros(0, nfe),
        pairwise_active = N_PAIR > 0,
        ethyl_acetate_soft = HAS_EA_SOFT,
        homotopy_log = homotopy_log,
        eps_final = isempty(homotopy_log) ? NaN : homotopy_log[end].eps,
    )

    return cStar, vStar, cdStar, lStar, alLStar, alUStar, auStar, apStar, aaaStar, aprStar, aeaStar, aatpmStar, hStar, diagStar
end
