# ======================================================================================
# SOLVER MULTIFASE CON OVERFLOW EN CRECIMIENTO
# ======================================================================================
# Arquitectura de 5 fases biologicas:
#   LAG             : max ATPM, sin crecimiento, sin turnover
#   GROWTH          : near-optimal mu + overflow (ATPM, isoamylol, isoamyl acetate)
#   GROWTH_NLIM     : growth bajo limitacion N, mayor sesgo overflow
#   STATIONARY      : max protein + ATPM, turnover activo, N_rec recycling
#   DECAY           : similar stationary con decay reforzado
# ======================================================================================

import math
import numpy as np
import pandas as pd
import cobra
import matplotlib.pyplot as plt

print("\n" + "=" * 100)
print("Solver multifase con overflow en crecimiento")
print("=" * 100)

# ── Aliases locales ──
_kl_fn_mp  = globals().get("_kinetic_limits_fn", globals().get("_kinetic_limits"))
_NV_IDS_MP = globals().get("N_VITAMIN_EX_IDS_LOC", globals().get("N_VITAMIN_EX_IDS", []))
_EPS_MP    = float(globals().get("EPS_LOC", globals().get("EPS", 1e-9)))
_VPROT_FF  = float(globals().get("VPROT_FLOOR_FRAC", 0.90))

# ── Parametros multifase (PSO los sobreescribe via globals) ──
ALPHA_GROWTH                = float(globals().get("ALPHA_GROWTH", 0.95))
W_ATPM_GROWTH               = float(globals().get("W_ATPM_GROWTH", 1.0))
W_ISOAMYLOL_GROWTH          = float(globals().get("W_ISOAMYLOL_GROWTH", 0.5))
W_ISOAMYL_ACETATE_GROWTH    = float(globals().get("W_ISOAMYL_ACETATE_GROWTH", 0.3))
ALPHA_GROWTH_NLIM           = float(globals().get("ALPHA_GROWTH_NLIM", 0.90))
LAG_END_H                   = float(globals().get("LAG_END_H", 4.0))
N_LIMITED_THRESHOLD          = float(globals().get("N_LIMITED_THRESHOLD", 0.02))
DECAY_START_H               = float(globals().get("DECAY_START_H", 60.0))

ISOAMYLOL_RID       = "r_1865"
ISOAMYL_ACETATE_RID = "r_1862"
_ATPM_GR_LB = float(globals().get("BAL_ATPM_GROWTH_LB", 0.7))
_ATPM_GR_UB = float(globals().get("BAL_ATPM_GROWTH_UB", 1000.0))

print(f"ALPHA_GROWTH={ALPHA_GROWTH}  W_ATPM={W_ATPM_GROWTH}  "
      f"W_ISO={W_ISOAMYLOL_GROWTH}  W_IAA={W_ISOAMYL_ACETATE_GROWTH}")
print(f"LAG_END_H={LAG_END_H}  N_LIM_THR={N_LIMITED_THRESHOLD}  DECAY_H={DECAY_START_H}")

# Collocation constants (same as cell 44)
_COLMAT = np.array([
    [0.19681547722366, -0.06553542585020,  0.02377097434822],
    [0.39442431473909,  0.29207341166523, -0.04154875212600],
    [0.37640306270047,  0.51248582618842,  0.11111111111111],
])
_COLMAT_INV = np.linalg.inv(_COLMAT)
_RADAU = np.array([0.15505102572168, 0.64494897427832, 1.0])


# ══════════════════════════════════════════════════════════════════════════════════════
# PHASE DETECTION
# ══════════════════════════════════════════════════════════════════════════════════════

def _infer_phase(cN_total, t_h):
    if t_h < LAG_END_H:
        return "lag"
    if cN_total <= N_TOTAL_DEPLETION_THRESHOLD:
        return "decay" if t_h >= DECAY_START_H else "stationary"
    if cN_total <= N_LIMITED_THRESHOLD:
        return "growth_N_limited"
    return "growth"


# ══════════════════════════════════════════════════════════════════════════════════════
# LP SOLVER MULTIFASE
# ══════════════════════════════════════════════════════════════════════════════════════

def _fba_at_state_colloc_multiphase(yk, aerobic_mode, t_h=0.0):
    yk = np.maximum(yk, 0.0)
    cX      = yk[IDX_X]
    cN_free = yk[IDX_NFREE]
    cG      = yk[IDX_G]
    cF      = yk[IDX_F]
    cE      = yk[IDX_E]
    cO2     = yk[IDX_O2]
    cProt   = yk[IDX_PROT]
    cCarb   = yk[IDX_CARB]
    cNrec   = yk[IDX_NREC]
    aa_state = {ak: float(yk[IDX_AA0 + j]) for j, ak in enumerate(aa_keys)}

    aa_n_sum = MW_N * sum(aa_state.values())
    cN_total = cN_free + aa_n_sum

    P_frac   = cProt / max(cX, _EPS_MP)
    C_frac   = cCarb / max(cX, _EPS_MP)
    full_gam = _compute_full_gam(P_frac, RNA_FRAC, C_frac,
                                  Pbase_global, Rbase_global, Cbase_global)
    lim_glu, lim_fru, lim_eth, lim_obj, lim_n_map = _kl_fn_mp(
        cX, cN_free, cG, cF, cE, cO2)

    phase = _infer_phase(cN_total, t_h)

    with model as mtmp:
        rb = {r.id: r for r in mtmp.reactions}

        for rid in _tracked_ids:
            if rid in rb:
                rb[rid].lower_bound = _base_lb_nrec[rid]
                rb[rid].upper_bound = _base_ub_nrec[rid]

        if GLU_ID in rb:
            rb[GLU_ID].lower_bound = max(rb[GLU_ID].lower_bound, -lim_glu)
            rb[GLU_ID].upper_bound = min(rb[GLU_ID].upper_bound,  0.0)
        if FRU_ID in rb:
            rb[FRU_ID].lower_bound = max(rb[FRU_ID].lower_bound, -lim_fru)
            rb[FRU_ID].upper_bound = min(rb[FRU_ID].upper_bound,  0.0)
        for rid in KINETIC_N_SOURCE_IDS:
            if rid in rb:
                Li = float(lim_n_map.get(rid, 0.0))
                rb[rid].lower_bound = max(rb[rid].lower_bound, -Li)
                rb[rid].upper_bound = min(rb[rid].upper_bound,  0.0)

        if aerobic_mode:
            lim_o2 = O2_VMAX_UPTAKE * (cO2 / (cO2 + O2_KO_G_L + _EPS_MP))
            lim_o2 = 0.0 if cO2 <= O2_DEPLETION_THRESHOLD else max(0.0, lim_o2)
            rb[O2_ID].lower_bound = -lim_o2
            rb[O2_ID].upper_bound =  0.0
        else:
            rb[O2_ID].lower_bound = 0.0
            rb[O2_ID].upper_bound = 0.0

        if APPLY_PRODUCT_CAPS:
            rb[ETH_ID].upper_bound = min(rb[ETH_ID].upper_bound, lim_eth)
            rb[OBJ_ID].upper_bound = min(rb[OBJ_ID].upper_bound, lim_obj)
        else:
            rb[ETH_ID].upper_bound = max(rb[ETH_ID].upper_bound, 1000.0)
            rb[OBJ_ID].upper_bound = max(rb[OBJ_ID].upper_bound, 1000.0)

        _apply_variable_gam(mtmp, full_gam)

        # ── LAG ──
        if phase == "lag":
            rb[OBJ_ID].lower_bound = 0.0
            rb[OBJ_ID].upper_bound = 0.0
            for rid in KINETIC_N_SOURCE_IDS:
                if rid in rb:
                    rb[rid].lower_bound = 0.0
                    rb[rid].upper_bound = 0.0
            for rid in _NV_IDS_MP:
                if rid in rb:
                    rb[rid].lower_bound = 0.0
                    rb[rid].upper_bound = 0.0
            rb[ATPM_ID].lower_bound = max(rb[ATPM_ID].lower_bound, ATPM_LB_NO_GROWTH)
            rb[ATPM_ID].upper_bound = ATPM_UB_NO_GROWTH
            mtmp.objective = ATPM_ID
            mtmp.objective_direction = "max"
            try:
                sol = cobra.flux_analysis.pfba(mtmp, fraction_of_optimum=1.0)
            except Exception:
                return None, "FAIL", phase
            mode = "LAG"

        # ── GROWTH / GROWTH_NLIM ──
        elif phase in ("growth", "growth_N_limited"):
            is_nlim = (phase == "growth_N_limited")
            _alpha = ALPHA_GROWTH_NLIM if is_nlim else ALPHA_GROWTH

            for ak, rid in found_precursors.items():
                if rid in rb:
                    q_i = max(0.0, K_AA_UPTAKE_GROWTH * aa_state[ak] / max(cX, _EPS_MP))
                    rb[rid].lower_bound = -q_i
                    rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)

            if ATPM_ID in rb:
                rb[ATPM_ID].lower_bound = max(rb[ATPM_ID].lower_bound, _ATPM_GR_LB)
                rb[ATPM_ID].upper_bound = max(rb[ATPM_ID].upper_bound, _ATPM_GR_UB)

            # Step A: max biomass
            mtmp.objective = OBJ_ID
            mtmp.objective_direction = "max"
            try:
                sol_mu = mtmp.optimize()
            except Exception:
                return None, "FAIL", phase
            if getattr(sol_mu, "status", "") != "optimal":
                return None, "NO_OPT", phase
            mu_star = max(0.0, float(sol_mu.fluxes.get(OBJ_ID, 0.0)))

            # Step B: epsilon-constraint
            rb[OBJ_ID].lower_bound = max(rb[OBJ_ID].lower_bound,
                                          float(_alpha) * mu_star)

            # Step C: max overflow
            _w_atpm = W_ATPM_GROWTH
            _w_iso  = W_ISOAMYLOL_GROWTH
            _w_iaa  = W_ISOAMYL_ACETATE_GROWTH
            if is_nlim:
                _w_iso *= 1.5
                _w_iaa *= 1.5

            mix_obj = {}
            if _w_atpm > 0 and ATPM_ID in rb:
                mix_obj[mtmp.reactions.get_by_id(ATPM_ID)] = float(_w_atpm)
            if _w_iso > 0 and ISOAMYLOL_RID in rb:
                mix_obj[mtmp.reactions.get_by_id(ISOAMYLOL_RID)] = float(_w_iso)
            if _w_iaa > 0 and ISOAMYL_ACETATE_RID in rb:
                mix_obj[mtmp.reactions.get_by_id(ISOAMYL_ACETATE_RID)] = float(_w_iaa)

            if mix_obj:
                mtmp.objective = mix_obj
                mtmp.objective_direction = "max"
                try:
                    sol = cobra.flux_analysis.pfba(mtmp, fraction_of_optimum=1.0)
                except Exception:
                    sol = sol_mu
            else:
                sol = sol_mu
            mode = "BIOMASS_NLIM" if is_nlim else "BIOMASS"

        # ── STATIONARY / DECAY ──
        else:
            for rid in KINETIC_N_SOURCE_IDS:
                if rid in rb:
                    rb[rid].lower_bound = max(rb[rid].lower_bound, 0.0)
                    rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)
            for rid in _NV_IDS_MP:
                if rid in rb:
                    rb[rid].lower_bound = max(rb[rid].lower_bound, 0.0)
                    rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)
            rb[OBJ_ID].lower_bound = 0.0
            rb[OBJ_ID].upper_bound = 0.0

            q_Nrec = max(0.0, K_NREC_UPTAKE * cNrec / max(cX, _EPS_MP))
            q_AA   = q_Nrec / MW_N
            for ak, rid in found_precursors.items():
                if rid in rb:
                    rb[rid].lower_bound = -max(0.0, q_AA * AA_ALPHA[ak])
                    rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)

            vprot_star = 0.0
            try:
                mtmp.objective = PROT_RXN_ID
                mtmp.objective_direction = "max"
                sol_p = mtmp.optimize()
                if getattr(sol_p, "status", "") == "optimal":
                    vprot_star = max(0.0, float(sol_p.fluxes.get(PROT_RXN_ID, 0.0)))
            except Exception:
                vprot_star = 0.0
            if vprot_star > 1e-10:
                rb[PROT_RXN_ID].lower_bound = max(rb[PROT_RXN_ID].lower_bound,
                                                   _VPROT_FF * vprot_star)
            rb[ATPM_ID].lower_bound = max(rb[ATPM_ID].lower_bound, ATPM_LB_NO_GROWTH)
            rb[ATPM_ID].upper_bound = ATPM_UB_NO_GROWTH
            mtmp.objective = ATPM_ID
            mtmp.objective_direction = "max"
            try:
                sol = cobra.flux_analysis.pfba(mtmp, fraction_of_optimum=1.0)
            except Exception:
                return None, "FAIL", phase
            mode = "DECAY" if phase == "decay" else "TURNOVER_ATPM"

        if getattr(sol, "status", "optimal") != "optimal":
            return None, "NO_OPT", phase

        fluxes = dict(
            v_obj=float(sol.fluxes.get(OBJ_ID, 0.0)),
            v_glu=float(sol.fluxes.get(GLU_ID, 0.0)),
            v_fru=float(sol.fluxes.get(FRU_ID, 0.0)),
            v_eth=float(sol.fluxes.get(ETH_ID, 0.0)),
            v_o2=float(sol.fluxes.get(O2_ID, 0.0)),
            v_atpm=float(sol.fluxes.get(ATPM_ID, 0.0)),
            v_prot=float(sol.fluxes.get(PROT_RXN_ID, 0.0)),
            full_gam=full_gam,
        )
        vn_eff = 0.0
        if mode in ("BIOMASS", "BIOMASS_NLIM"):
            for rid in KINETIC_N_SOURCE_IDS:
                vi = float(sol.fluxes.get(rid, 0.0))
                vn_eff += max(0.0, -vi) * N_atoms_map.get(rid, 1.0) * MW_N
        fluxes["vn_eff"] = vn_eff

        for ak in aa_keys:
            fluxes[f"v_aa_{ak}"] = float(sol.fluxes.get(found_precursors[ak], 0.0))
        for ak_a in aroma_keys:
            fluxes[f"v_aroma_{ak_a}"] = float(sol.fluxes.get(found_aromas[ak_a], 0.0))
        for _rid in globals().get("EXTRA_AROMA_RIDS", []):
            fluxes[f"v_aroma_{_rid}"] = float(sol.fluxes.get(_rid, 0.0))

    return fluxes, mode, phase


# ══════════════════════════════════════════════════════════════════════════════════════
# RHS CON TURNOVER CONDICIONAL A FASE
# ══════════════════════════════════════════════════════════════════════════════════════

def _compute_rhs_colloc_multiphase(yk, fluxes, mode, phase):
    cX    = max(0.0, yk[IDX_X])
    cProt = max(0.0, yk[IDX_PROT])
    cCarb = max(0.0, yk[IDX_CARB])
    cNrec = max(0.0, yk[IDX_NREC])
    Xa    = XA_FRACTION * cX

    rhs = np.zeros(n_state)
    is_growth = mode in ("BIOMASS", "BIOMASS_NLIM")

    rhs[IDX_X]     = max(0.0, fluxes["v_obj"]) * cX if is_growth else 0.0
    rhs[IDX_NFREE] = -fluxes["vn_eff"] * cX if is_growth else 0.0
    rhs[IDX_G]     = -MW_GLU * max(0.0, -fluxes["v_glu"]) * cX
    rhs[IDX_F]     = -MW_FRU * max(0.0, -fluxes["v_fru"]) * cX
    rhs[IDX_E]     =  MW_ETH * max(0.0,  fluxes["v_eth"]) * cX
    rhs[IDX_O2]    = -MW_O2  * max(0.0, -fluxes["v_o2"])  * cX

    # Turnover SOLO en stationary/decay
    turnover_active = phase in ("stationary", "decay")
    lam = TURNOVER_LAMBDA if turnover_active else 0.0

    rhs[IDX_PROT] = (
        max(0.0, fluxes["v_obj"]) * cX * PROT_CONTENT_0
        - cProt * K_DEATH
        + Xa * max(0.0, fluxes["v_prot"])
        - lam * cProt
    )
    rhs[IDX_CARB] = (
        max(0.0, fluxes["v_obj"]) * cX * CARB_CONTENT_0
        - cCarb * K_DEATH
    )

    nrec_inflow = Y_N_FROM_PROT * lam * cProt
    if mode in ("TURNOVER_ATPM", "DECAY"):
        total_aa_mmol = sum(max(0.0, -fluxes.get(f"v_aa_{ak}", 0.0)) for ak in aa_keys)
        nrec_outflow = total_aa_mmol * MW_N * cX
    else:
        nrec_outflow = 0.0
    rhs[IDX_NREC] = nrec_inflow - nrec_outflow

    for j_aa, ak in enumerate(aa_keys):
        v_aa = fluxes.get(f"v_aa_{ak}", 0.0)
        rhs[IDX_AA0 + j_aa] = cX * v_aa if is_growth else 0.0

    for j_ar, ak_a in enumerate(aroma_keys):
        mw_a = AROMA_MW.get(ak_a, 0.1)
        rhs[IDX_AROMA0 + j_ar] = mw_a * max(0.0, fluxes.get(f"v_aroma_{ak_a}", 0.0)) * cX

    return rhs


# ══════════════════════════════════════════════════════════════════════════════════════
# FUNCION PRINCIPAL — dFBA SOA multifase con Radau IIA
# ══════════════════════════════════════════════════════════════════════════════════════

def _run_dfba_colloc_multiphase(
    *,
    scenario_name="Multifase + overflow",
    aerobic_mode=False,
    o2_init_g_l=0.0,
    nfe=18,
    ncp=3,
    t_end=72.0,
    n_iter=2,
):
    assert ncp == 3
    h = t_end / nfe
    hm = np.full(nfe, h)
    nc = n_state

    t_fe_start = np.concatenate([[0.0], np.cumsum(hm[:-1])])
    t_loc = np.array([[t_fe_start[i] + _RADAU[j] * h
                       for j in range(ncp)] for i in range(nfe)])

    C    = np.zeros((nc, nfe, ncp))
    CDOT = np.zeros((nc, nfe, ncp))
    V_flux  = [[None] * ncp for _ in range(nfe)]
    V_mode  = [[""]   * ncp for _ in range(nfe)]
    V_phase = [[""]   * ncp for _ in range(nfe)]

    N0_total   = float(globals().get("N0", 0.14))
    N0_ammonia = N0_total * N_AMMONIA_FRACTION
    N0_from_aa = max(0.0, N0_total - N0_ammonia)
    aa0_each   = (N0_from_aa / MW_N) / max(1, n_aa)

    c0 = np.zeros(nc)
    c0[IDX_X]     = float(globals().get("X0_init", 0.5))
    c0[IDX_NFREE] = N0_ammonia
    c0[IDX_G]     = float(globals().get("G0_init", 110.0))
    c0[IDX_F]     = float(globals().get("F0_init", 110.0))
    c0[IDX_E]     = float(globals().get("E0_init", 0.0))
    c0[IDX_O2]    = o2_init_g_l
    c0[IDX_PROT]  = c0[IDX_X] * PROT_CONTENT_0
    c0[IDX_CARB]  = c0[IDX_X] * CARB_CONTENT_0
    c0[IDX_NREC]  = N_REC_INIT
    for j in range(n_aa):
        c0[IDX_AA0 + j] = aa0_each

    c_prev = c0.copy()
    C_endpoints = np.zeros((nfe + 1, nc))
    C_endpoints[0] = c0

    print(f"\n[{scenario_name}]  nfe={nfe}, ncp={ncp}, h={h:.3f} h, n_iter={n_iter}")

    for i in range(nfe):
        for j in range(ncp):
            C[:, i, j] = np.maximum(c_prev, 0.0)

        for _it in range(n_iter):
            for j in range(ncp):
                _t_cp = float(t_loc[i, j])

                _tp = globals().get("T_PROFILE_DF", None)
                if _tp is not None and len(_tp):
                    _T_now = float(np.interp(_t_cp,
                                             _tp["t_h"].to_numpy(),
                                             _tp["T_K"].to_numpy()))
                else:
                    _T_now = float(dynamic_temperature(_t_cp))
                globals()["T_val"] = _T_now

                result = _fba_at_state_colloc_multiphase(
                    C[:, i, j], aerobic_mode, t_h=_t_cp)

                if result[0] is None:
                    fluxes = V_flux[i][j]
                    mode   = "FAIL"
                    phase  = result[2] if len(result) > 2 else "unknown"
                    if fluxes is None:
                        zero_f = {k: 0.0 for k in [
                            "v_obj","v_glu","v_fru","v_eth",
                            "v_o2","v_atpm","v_prot","vn_eff","full_gam"]}
                        for ak in aa_keys:
                            zero_f[f"v_aa_{ak}"] = 0.0
                        for ak_a in aroma_keys:
                            zero_f[f"v_aroma_{ak_a}"] = 0.0
                        for _rid_extra in globals().get("EXTRA_AROMA_RIDS", []):
                            zero_f[f"v_aroma_{_rid_extra}"] = 0.0
                        fluxes = zero_f
                else:
                    fluxes, mode, phase = result

                V_flux[i][j]  = fluxes
                V_mode[i][j]  = mode
                V_phase[i][j] = phase

                CDOT[:, i, j] = _compute_rhs_colloc_multiphase(
                    C[:, i, j], fluxes, mode, phase)

            C_new = c_prev[:, None] + h * (CDOT[:, i, :] @ _COLMAT.T)
            C[:, i, :] = np.maximum(C_new, 0.0)

        c_prev = C[:, i, ncp - 1].copy()

        _t_start_fe = h * i
        _t_end_fe   = h * (i + 1)
        _pulses = globals().get("PULSE_SCHEDULE", None) or []
        for _p in _pulses:
            _tp_pulse = float(_p["t_h"])
            if _t_start_fe < _tp_pulse <= _t_end_fe:
                c_prev[IDX_NFREE] += float(_p.get("amount_gN_L", 0.0))

        C_endpoints[i + 1] = c_prev

    # ODE gap
    max_ode_gap = 0.0
    for i in range(nfe):
        c_p = C_endpoints[i]
        for s in range(nc):
            rhs_j = np.array([(C[s, i, j] - c_p[s]) / h for j in range(ncp)])
            cdot_coll = _COLMAT_INV @ rhs_j
            gap = float(np.max(np.abs(cdot_coll - CDOT[s, i, :])))
            if gap > max_ode_gap:
                max_ode_gap = gap

    # ── Unpack ──
    t_all = np.concatenate([[0.0], t_loc.ravel()])
    n_pts = 1 + nfe * ncp
    Y_all = np.zeros((n_pts, nc))
    Y_all[0] = c0
    for i in range(nfe):
        for j in range(ncp):
            Y_all[1 + i * ncp + j] = C[:, i, j]

    t_colloc = t_loc.ravel()

    def _fs(key):
        out = np.zeros(nfe * ncp)
        for i in range(nfe):
            for j in range(ncp):
                f = V_flux[i][j]
                out[i * ncp + j] = f.get(key, 0.0) if f else 0.0
        return out

    mu_c    = _fs("v_obj")
    vglu_c  = np.maximum(0.0, -_fs("v_glu"))
    vfru_c  = np.maximum(0.0, -_fs("v_fru"))
    veth_c  = np.maximum(0.0,  _fs("v_eth"))
    vo2_c   = np.maximum(0.0, -_fs("v_o2"))
    vn_c    = np.maximum(0.0,  _fs("vn_eff"))
    vatpm_c = np.maximum(0.0,  _fs("v_atpm"))
    vprot_c = np.maximum(0.0,  _fs("v_prot"))
    gam_c   = _fs("full_gam")
    mode_c  = [V_mode[i][j]  for i in range(nfe) for j in range(ncp)]
    phase_c = [V_phase[i][j] for i in range(nfe) for j in range(ncp)]

    aa_uptake_c  = {ak: np.maximum(0.0, -_fs(f"v_aa_{ak}"))    for ak in aa_keys}
    aroma_flux_c = {ak: np.maximum(0.0,  _fs(f"v_aroma_{ak}")) for ak in aroma_keys}

    stat_mask = np.array([m in ("TURNOVER_ATPM", "DECAY") for m in mode_c])
    t_stat = float(t_colloc[np.argmax(stat_mask)]) if stat_mask.any() else float("nan")

    X_all     = Y_all[:, IDX_X]
    Nfree_all = Y_all[:, IDX_NFREE]
    Prot_all  = Y_all[:, IDX_PROT]
    Carb_all  = Y_all[:, IDX_CARB]
    Nrec_all  = Y_all[:, IDX_NREC]

    aa_conc_all = {ak: Y_all[:, IDX_AA0 + j]     for j,   ak  in enumerate(aa_keys)}
    aroma_all   = {ak: Y_all[:, IDX_AROMA0 + j_a] for j_a, ak in enumerate(aroma_keys)}

    aa_n_all      = (MW_N * np.sum(np.column_stack([aa_conc_all[ak] for ak in aa_keys]), axis=1)
                     if aa_keys else np.zeros(n_pts))
    N_in_prot_all = Y_N_FROM_PROT * Prot_all
    N_total_all   = Nfree_all + aa_n_all + N_in_prot_all + Nrec_all
    N_bal_err     = abs(float(N_total_all[-1]) - float(N_total_all[0]))

    P_frac_all = Prot_all / np.maximum(X_all, _EPS_MP)
    C_frac_all = Carb_all / np.maximum(X_all, _EPS_MP)
    gam_full_all = np.array([
        _compute_full_gam(float(P_frac_all[ii]), RNA_FRAC, float(C_frac_all[ii]),
                          Pbase_global, Rbase_global, Cbase_global)
        for ii in range(n_pts)
    ])

    # Phase counts
    phase_counts = {}
    for p in phase_c:
        phase_counts[p] = phase_counts.get(p, 0) + 1

    print(f"  Phase counts: {phase_counts}")
    print(f"  t_stationary ~ {t_stat:.1f} h")
    print(f"  X_final={X_all[-1]:.4f}  E_final={Y_all[-1, IDX_E]:.4f}")
    print(f"  N_bal: {N_total_all[0]:.6f} -> {N_total_all[-1]:.6f} "
          f"(err={N_bal_err:.2e}  {'OK' if N_bal_err < 0.01 else 'CHECK!'})")
    print(f"  ODE gap max = {max_ode_gap:.3e}")

    # Isoamylol / isoamyl acetate diagnostics
    _iso_key = None
    for _ak, _rid in found_aromas.items():
        if _rid == ISOAMYLOL_RID:
            _iso_key = _ak
            break
    if _iso_key and _iso_key in aroma_all:
        print(f"  Isoamylol ({_iso_key}) final = {1000*aroma_all[_iso_key][-1]:.4f} mg/L")
    _v_iaa = np.maximum(0.0, _fs(f"v_aroma_{ISOAMYL_ACETATE_RID}"))
    _iaa_max = float(np.max(_v_iaa)) if len(_v_iaa) else 0.0
    print(f"  Isoamyl acetate flux max = {_iaa_max:.6g} mmol/gDW/h")

    # ── DataFrames ──
    states_df = pd.DataFrame({
        "t_h":                 t_all,
        "X_gDW_L":             X_all,
        "N_gN_L":              Nfree_all + aa_n_all,
        "N_free_gN_L":         Nfree_all,
        "N_from_AA_gN_L":      aa_n_all,
        "N_rec_gN_L":          Nrec_all,
        "N_in_Prot_gN_L":      N_in_prot_all,
        "N_total_system_gN_L": N_total_all,
        "G_g_L":               Y_all[:, IDX_G],
        "F_g_L":               Y_all[:, IDX_F],
        "E_g_L":               Y_all[:, IDX_E],
        "O2_g_L":              Y_all[:, IDX_O2],
        "Prot_g_L":            Prot_all,
        "Carb_g_L":            Carb_all,
        "P_frac_gProt_gDW":    P_frac_all,
        "C_frac_gCarb_gDW":    C_frac_all,
        "GAM_mmol_gDW":        gam_full_all,
        "fe_idx":              np.concatenate([[-1], np.repeat(np.arange(nfe), ncp)]),
        "cp_idx":              np.concatenate([[-1], np.tile(np.arange(ncp),   nfe)]),
    })
    for ak in aa_keys:
        states_df[f"AA_{ak}_mmol_L"] = aa_conc_all[ak]
        states_df[f"AA_{ak}_g_L"]    = aa_conc_all[ak] * AA_MW.get(ak, 0.13)
    for ak in aroma_keys:
        states_df[f"{ak}_g_L"]  = aroma_all[ak]
        states_df[f"{ak}_mg_L"] = 1000.0 * aroma_all[ak]

    flux_df = pd.DataFrame({
        "t_h":                  t_colloc,
        "fe_idx":               np.repeat(np.arange(nfe), ncp),
        "cp_idx":               np.tile(np.arange(ncp), nfe),
        "mode":                 mode_c,
        "phase":                phase_c,
        "mu_h_inv":             mu_c,
        "v_glu":                vglu_c,
        "v_fru":                vfru_c,
        "v_eth":                veth_c,
        "v_o2":                 vo2_c,
        "v_N_uptake_gN_gDW_h":  vn_c,
        "v_ATPM":               vatpm_c,
        "v_Prot":               vprot_c,
        "GAM":                  gam_c,
    })
    for ak in aa_keys:
        flux_df[f"v_uptake_{ak}"] = aa_uptake_c[ak]
    for ak in aroma_keys:
        flux_df[f"v_{ak}"] = aroma_flux_c[ak]
    for _rid in globals().get("EXTRA_AROMA_RIDS", []):
        flux_df[f"v_aroma_{_rid}"] = np.maximum(0.0, _fs(f"v_aroma_{_rid}"))

    display(flux_df.head(9))
    display(flux_df.tail(9))

    # ── Plots ──
    _phase_colors = {"lag": "gray", "growth": "green", "growth_N_limited": "olive",
                     "stationary": "orange", "decay": "red"}

    state_series = [
        ("X (gDW/L)",      X_all,              np.interp(t_all, t_colloc, mu_c),    "mu (1/h)"),
        ("N_ext (gN/L)",   Nfree_all+aa_n_all, np.interp(t_all, t_colloc, vn_c),    "v_N"),
        ("Glucosa (g/L)",  Y_all[:, IDX_G],    np.interp(t_all, t_colloc, vglu_c),  "v_glu"),
        ("Fructosa (g/L)", Y_all[:, IDX_F],    np.interp(t_all, t_colloc, vfru_c),  "v_fru"),
        ("Etanol (g/L)",   Y_all[:, IDX_E],    np.interp(t_all, t_colloc, veth_c),  "v_eth"),
        ("Prot (g/L)",     Prot_all,           np.interp(t_all, t_colloc, vprot_c), "v_Prot"),
        ("N_rec (gN/L)",   Nrec_all,           np.zeros(n_pts),                      ""),
        ("N total (gN/L)", N_total_all,        np.zeros(n_pts),                      "balance"),
    ]
    ncols_s = 4
    nrows_s = math.ceil(len(state_series) / ncols_s)
    fig, ax = plt.subplots(nrows_s, ncols_s,
                            figsize=(5.2*ncols_s, 3.5*nrows_s), sharex=True, dpi=120)
    ax = np.atleast_1d(ax).ravel()
    for i_p, (ttl, yy_, vv_, vlbl) in enumerate(state_series):
        ax[i_p].plot(t_all, yy_, lw=2, color="tab:blue")
        # Phase background
        for k_cp in range(len(t_colloc)):
            pc = _phase_colors.get(phase_c[k_cp], "white")
            t0 = t_colloc[k_cp - 1] if k_cp > 0 else 0.0
            t1 = t_colloc[k_cp]
            ax[i_p].axvspan(t0, t1, alpha=0.08, color=pc)
        if np.isfinite(t_stat):
            ax[i_p].axvline(t_stat, ls="--", lw=1, color="navy", alpha=0.5)
        ax[i_p].set_title(ttl, fontsize=9)
        ax[i_p].grid(alpha=0.3)
        ax2 = ax[i_p].twinx()
        vv_cp = vv_ if len(vv_) == len(t_colloc) else np.interp(t_colloc, t_all, vv_)
        ax2.step(t_colloc, vv_cp, where="post", lw=1.4, ls="--",
                 color="tab:orange", alpha=0.85)
        ax2.set_ylabel(vlbl, color="tab:orange", fontsize=7)
    for j_p in range(i_p + 1, len(ax)):
        ax[j_p].axis("off")
    plt.suptitle(f"Multifase [{nfe} FE, {ncp} CP] — {scenario_name}", fontsize=11, y=1.01)
    plt.tight_layout()
    plt.show()

    # Aromas
    if aroma_keys:
        ncols_a = min(len(aroma_keys), 5)
        nrows_a = math.ceil(len(aroma_keys) / ncols_a)
        fig_a, ax_a = plt.subplots(nrows_a, ncols_a,
                                    figsize=(4.5*ncols_a, 3.5*nrows_a), dpi=120, squeeze=False)
        ax_a = ax_a.ravel()
        for i_a, ak_a in enumerate(aroma_keys):
            ax_a[i_a].plot(t_all, 1000.0*aroma_all[ak_a], lw=2.5, color="tab:brown")
            ax_a[i_a].set_title(f"{ak_a} (mg/L)", fontsize=9)
            ax_a[i_a].set_xlabel("t (h)")
            ax_a[i_a].grid(alpha=0.3)
        for j_a in range(i_a + 1, len(ax_a)):
            ax_a[j_a].axis("off")
        plt.suptitle(f"Aromas – {scenario_name}", fontsize=11)
        plt.tight_layout()
        plt.show()

    sfx = "aerobic_mp" if aerobic_mode else "anaerobic_mp"
    out_s = OUT_DIR / f"dfba_colloc_states_{sfx}.csv"
    out_f = OUT_DIR / f"dfba_colloc_flux_{sfx}.csv"
    states_df.to_csv(out_s, index=False)
    flux_df.to_csv(out_f, index=False)
    print(f"\nArchivos: {out_s}\n  {out_f}")

    return dict(
        states_df=states_df, flux_df=flux_df,
        C=C, CDOT=CDOT, C_endpoints=C_endpoints,
        V_flux=V_flux, V_mode=V_mode, V_phase=V_phase,
        t_loc=t_loc, hm=hm, t_stat=t_stat,
        nfe=nfe, ncp=ncp, max_ode_gap=max_ode_gap,
    )


# ══════════════════════════════════════════════════════════════════════════════════════
# DEMO
# ══════════════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("Demo: Multifase + overflow  (nfe=18, ncp=3, n_iter=2)")
print("=" * 80)

res_mp = _run_dfba_colloc_multiphase(
    scenario_name="Multifase + overflow",
    aerobic_mode=False,
    o2_init_g_l=0.0,
    nfe=18, ncp=3, t_end=72.0, n_iter=2,
)

# Update global references for downstream cells
states_df_phase = res_mp["states_df"].copy()
flux_df_phase = res_mp["flux_df"].copy()
print("[INFO] states_df_phase y flux_df_phase actualizados con trayectoria multifase.")
