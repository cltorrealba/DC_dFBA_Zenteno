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
import MathOptInterface as MOI

_softplus(x, ϵ) = 0.5 * (x + sqrt(x * x + ϵ * ϵ))
_sigmoid(x) = 1.0 / (1.0 + exp(-x))

function pFBA_KKT_flux_Zenteno_vargam_simultaneous(
    c0;
    eps_flux::Float64 = 0.0,
    apply_product_caps::Bool = false,
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

    # ---------- Radau IIA ----------
    colmat = [0.19681547722366  -0.06553542585020   0.02377097434822;
              0.39442431473909   0.29207341166523  -0.04154875212600;
              0.37640306270047   0.51248582618842   0.11111111111111]
    radau  = [0.15505102572168, 0.64494897427832, 1.0]

    # ---------- Model ----------
    m = Model(optimizer_with_attributes(
        Ipopt.Optimizer,
        "linear_solver" => IPOPT_LINEAR_SOLVER,
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
        "warm_start_init_point" => "yes"
    ))

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
    for i in 1:nfe
        set_start_value(hv[i], hm[i])
        for j in 1:ncp
            for s in 1:nc
                set_start_value(c[s, i, j], c0[s] / cs[s])
            end
        end
        for rx in 1:nv
            set_start_value(v[rx, i], 0.0)
        end
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

        Sc[mc=1:nm, i=1:nfe],
            sum(S[mc,k] * v[k,i] * vs[k] for k in 1:nv) == 0

        v_UB[mc=1:nv, i=1:nfe],  v[mc,i] * vs[mc] - vub[mc] <= 0
        v_LB[mc=1:nv, i=1:nfe], -v[mc,i] * vs[mc] + vlb[mc] <= 0

        MFE1, sum(hv[i] for i in 1:nfe) == th
        MFE3[i=1:nfe], hv[i] >= (1.0 - var_h) * hm[1]
        MFE4[i=1:nfe], hv[i] <= (1.0 + var_h) * hm[1]

        flux_smooth_pos[mc=1:nv, i=2:nfe],  v[mc,i] - v[mc,i-1] <= eps_flux
        flux_smooth_neg[mc=1:nv, i=2:nfe],  v[mc,i-1] - v[mc,i] <= eps_flux
    end)

    # ---------- stationarity (w.r.t. fluxes) ----------
    @NLconstraint(m, lagr_stationarity[mc=1:nv, i=1:nfe],
        d_phase[mc,i] +
        Q_REG * (v[mc,i] * vs[mc]) +
        alpha_L[mc,i] + alpha_U[mc,i] +
        sum(SELECT_UPTAKE[mc,k] * alpha_upt[k,i] for k in 1:N_UP) +
        sum(SELECT_PRODUCT[mc,k] * alpha_prod[k,i] for k in 1:N_PROD) +
        sum(SELECT_AAUPTAKE[mc,a] * alpha_aa[a,i] for a in 1:N_AA) +
        sum((PAIR_ALCOHOL_SELECT[mc,p] * PAIRWISE_PHI[p] - PAIR_ESTER_SELECT[mc,p]) * alpha_pair[p,i]
            for p in 1:N_PAIR) +
        (EA_ALCOHOL_SELECT[mc] * PHI_ETHYL_ACETATE_STATIC - EA_ESTER_SELECT[mc]) * alpha_ea[i] +
        (OBJ_SELECT[mc] * GAM_EXTRA_FE[i] - ATPM_SELECT[mc]) * alpha_atpm_floor[i] +
        sum(S[k,mc] * lambda_[k,i] for k in 1:nm) == 0
    )

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
