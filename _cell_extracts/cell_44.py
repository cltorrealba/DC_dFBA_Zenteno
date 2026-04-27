import math
import re
import numpy as np
import pandas as pd
import cobra
import matplotlib.pyplot as plt

# ======================================================================================
# dFBA SOA + pFBA con colocaciones ortogonales Radau IIA (ncp=3)
# Emula la estructura del modelo simultaneo de Notebook 2, pero resuelto
# de forma secuencial (SOA): ncp LPs por elemento finito, con pasos de
# estado definidos por la formula de colocacion de Radau IIA.
# ======================================================================================

print("\n" + "=" * 100)
print("dFBA SOA + pFBA con GAM variable + N_rec + Colocaciones ortogonales Radau IIA")
print("=" * 100)

# Verificar dependencias
_required_colloc = [
    "model", "rxn_by_id", "OBJ_ID", "GLU_ID", "FRU_ID", "ETH_ID", "O2_ID", "ATPM_ID",
    "KINETIC_N_SOURCE_IDS", "N_VITAMIN_EX_IDS", "N_atoms_map",
    "ATPM_LB_NO_GROWTH", "ATPM_UB_NO_GROWTH",
    "_kinetic_limits", "DFBA_DT",
    "N_DEPLETION_CLOSE_THRESHOLD", "O2_DEPLETION_THRESHOLD",
    "O2_VMAX_UPTAKE", "O2_KO_G_L",
    "MW_GLU", "MW_FRU", "MW_ETH", "MW_N", "MW_O2",
    "OUT_DIR", "APPLY_PRODUCT_CAPS", "EPS",
    "_compute_full_gam", "_apply_variable_gam",
    "Pbase_global", "Cbase_global", "Rbase_global", "RNA_FRAC",
    "found_precursors", "found_aromas", "PROT_RXN_ID",
    "AA_MW", "AROMA_MW",
    "PROT_CONTENT_0", "CARB_CONTENT_0", "K_DEATH", "TURNOVER_LAMBDA", "XA_FRACTION",
    "N_AMMONIA_FRACTION", "N_TOTAL_DEPLETION_THRESHOLD",
    "K_AA_UPTAKE_GROWTH", "BIOMASS_LOCK_FRAC", "AA_UPTAKE_WEIGHT",
    "Y_N_FROM_PROT", "K_NREC_UPTAKE", "N_REC_INIT",
    "AA_ALPHA", "aa_keys", "aroma_keys", "n_aa", "n_aroma",
    "_tracked_ids", "_base_lb_nrec", "_base_ub_nrec",
]
_missing_colloc = [s for s in _required_colloc if s not in globals()]
if _missing_colloc:
    raise RuntimeError(f"Faltan simbolos requeridos (ejecutar celda 44 primero): {_missing_colloc}")

_kinetic_limits_fn_c = globals().get("_kinetic_limits_fn", globals().get("_kinetic_limits"))
N_VITAMIN_EX_IDS_C   = globals().get("N_VITAMIN_EX_IDS_LOC", globals().get("N_VITAMIN_EX_IDS", []))
EPS_C                = float(globals().get("EPS_LOC", globals().get("EPS", 1e-9)))

# Rebuild state indexing locally to avoid stale global indices from previous runs.
IDX_X      = 0
IDX_NFREE  = 1
IDX_G      = 2
IDX_F      = 3
IDX_E      = 4
IDX_O2     = 5
IDX_PROT   = 6
IDX_CARB   = 7
IDX_NREC   = 8
IDX_AA0    = 9
IDX_AROMA0 = 9 + n_aa
n_state    = 9 + n_aa + n_aroma

# ══════════════════════════════════════════════════════════════════════════════════════
# CONSTANTES RADAU IIA (ncp=3) — identicas a Notebook 2
# ══════════════════════════════════════════════════════════════════════════════════════
#
# Butcher tableau A de Radau IIA, 3 etapas:
#   c[s,i,j] = c_prev[s] + h * sum_k  COLMAT[j,k] * cdot[s,i,k]
#   cdot      = COLMAT_INV @ [(C - c_prev) / h]   (CDOT colloc-consistente)
#
# Raices de Radau sobre [0,1]:  tau = [0.1550..., 0.6449..., 1.0]
# El ultimo punto (tau=1) coincide con el extremo derecho del FE,
# garantizando que el esquema es stiffly accurate y A-estable.

COLMAT_RADAU = np.array([
    [0.19681547722366, -0.06553542585020,  0.02377097434822],
    [0.39442431473909,  0.29207341166523, -0.04154875212600],
    [0.37640306270047,  0.51248582618842,  0.11111111111111],
])
COLMAT_RADAU_INV = np.linalg.inv(COLMAT_RADAU)
RADAU_ROOTS      = np.array([0.15505102572168, 0.64494897427832, 1.0])

print(f"Radau IIA  cond(COLMAT) = {np.linalg.cond(COLMAT_RADAU):.3f}")
print(f"Puntos de colocacion en FE: tau = {RADAU_ROOTS}")

# ══════════════════════════════════════════════════════════════════════════════════════
# HELPER 1 — LP / pFBA en un punto de colocacion
# ══════════════════════════════════════════════════════════════════════════════════════

def _fba_at_state_colloc(yk, aerobic_mode: bool):
    """
    Resuelve LP (FBA / pFBA) en el estado fisico yk.
    Misma logica de conmutacion de fase que celda 44 (_run_dfba_nrec):
      BIOMASS       -> maximiza crecimiento + AA uptake (pFBA secundario)
      TURNOVER_ATPM -> maximiza ATPM con pool N_rec

    Retorna:
        fluxes : dict con todos los flujos necesarios para el RHS
        mode   : 'BIOMASS' | 'TURNOVER_ATPM' | 'FAIL' | 'NO_OPT'
    """
    yk = np.maximum(yk, 0.0)
    cX       = yk[IDX_X]
    cN_free  = yk[IDX_NFREE]
    cG       = yk[IDX_G]
    cF       = yk[IDX_F]
    cE       = yk[IDX_E]
    cO2      = yk[IDX_O2]
    cProt    = yk[IDX_PROT]
    cCarb    = yk[IDX_CARB]
    cNrec    = yk[IDX_NREC]
    aa_state = {ak: float(yk[IDX_AA0 + j]) for j, ak in enumerate(aa_keys)}

    aa_n_sum = MW_N * sum(aa_state.values())
    cN_total = cN_free + aa_n_sum

    P_frac   = cProt / max(cX, EPS_C)
    C_frac   = cCarb / max(cX, EPS_C)
    full_gam = _compute_full_gam(P_frac, RNA_FRAC, C_frac,
                                  Pbase_global, Rbase_global, Cbase_global)

    lim_glu, lim_fru, lim_eth, lim_obj, lim_n_map = _kinetic_limits_fn_c(
        cX, cN_free, cG, cF, cE, cO2)

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
            lim_o2 = O2_VMAX_UPTAKE * (cO2 / (cO2 + O2_KO_G_L + EPS_C))
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

        # CONMUTACION DE FASE
        if cN_total > N_TOTAL_DEPLETION_THRESHOLD:
            # FASE CRECIMIENTO
            for ak, rid in found_precursors.items():
                if rid in rb:
                    q_i = max(0.0, K_AA_UPTAKE_GROWTH * aa_state[ak] / max(cX, EPS_C))
                    rb[rid].lower_bound = -q_i
                    rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)

            mtmp.objective           = OBJ_ID
            mtmp.objective_direction = "max"
            try:
                sol_mu = mtmp.optimize()
            except Exception:
                return None, "FAIL"
            if getattr(sol_mu, "status", "") != "optimal":
                return None, "NO_OPT"

            mu_star = max(0.0, float(sol_mu.fluxes.get(OBJ_ID, 0.0)))
            rb[OBJ_ID].lower_bound = max(rb[OBJ_ID].lower_bound,
                                          BIOMASS_LOCK_FRAC * mu_star)

            uptake_obj = {mtmp.reactions.get_by_id(rid): -AA_UPTAKE_WEIGHT
                          for rid in found_precursors.values() if rid in rb}
            mtmp.objective           = uptake_obj
            mtmp.objective_direction = "max"
            try:
                sol = cobra.flux_analysis.pfba(mtmp, fraction_of_optimum=1.0)
            except Exception:
                return None, "FAIL"
            mode = "BIOMASS"

        else:
            # FASE ESTACIONARIA
            for rid in KINETIC_N_SOURCE_IDS:
                if rid in rb:
                    rb[rid].lower_bound = max(rb[rid].lower_bound, 0.0)
                    rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)
            for rid in N_VITAMIN_EX_IDS_C:
                if rid in rb:
                    rb[rid].lower_bound = max(rb[rid].lower_bound, 0.0)
                    rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)

            rb[OBJ_ID].lower_bound = 0.0
            rb[OBJ_ID].upper_bound = 0.0

            q_Nrec = max(0.0, K_NREC_UPTAKE * cNrec / max(cX, EPS_C))
            q_AA   = q_Nrec / MW_N
            for ak, rid in found_precursors.items():
                if rid in rb:
                    rb[rid].lower_bound = -max(0.0, q_AA * AA_ALPHA[ak])
                    rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)

            vprot_star = 0.0
            try:
                mtmp.objective           = PROT_RXN_ID
                mtmp.objective_direction = "max"
                sol_p = mtmp.optimize()
                if getattr(sol_p, "status", "") == "optimal":
                    vprot_star = max(0.0, float(sol_p.fluxes.get(PROT_RXN_ID, 0.0)))
            except Exception:
                vprot_star = 0.0

            if vprot_star > 1e-10:
                rb[PROT_RXN_ID].lower_bound = max(rb[PROT_RXN_ID].lower_bound,
                                                   VPROT_FLOOR_FRAC * vprot_star)

            rb[ATPM_ID].lower_bound = max(rb[ATPM_ID].lower_bound, ATPM_LB_NO_GROWTH)
            rb[ATPM_ID].upper_bound = ATPM_UB_NO_GROWTH

            mtmp.objective           = ATPM_ID
            mtmp.objective_direction = "max"
            try:
                sol = cobra.flux_analysis.pfba(mtmp, fraction_of_optimum=1.0)
            except Exception:
                return None, "FAIL"
            mode = "TURNOVER_ATPM"

        if getattr(sol, "status", "optimal") != "optimal":
            return None, "NO_OPT"

        fluxes = dict(
            v_obj    = float(sol.fluxes.get(OBJ_ID,      0.0)),
            v_glu    = float(sol.fluxes.get(GLU_ID,      0.0)),
            v_fru    = float(sol.fluxes.get(FRU_ID,      0.0)),
            v_eth    = float(sol.fluxes.get(ETH_ID,      0.0)),
            v_o2     = float(sol.fluxes.get(O2_ID,       0.0)),
            v_atpm   = float(sol.fluxes.get(ATPM_ID,     0.0)),
            v_prot   = float(sol.fluxes.get(PROT_RXN_ID, 0.0)),
            full_gam = full_gam,
        )
        vn_eff = 0.0
        if mode == "BIOMASS":
            for rid in KINETIC_N_SOURCE_IDS:
                vi = float(sol.fluxes.get(rid, 0.0))
                vn_eff += max(0.0, -vi) * N_atoms_map.get(rid, 1.0) * MW_N
        fluxes["vn_eff"] = vn_eff

        for j_aa, ak in enumerate(aa_keys):
            fluxes[f"v_aa_{ak}"] = float(sol.fluxes.get(found_precursors[ak], 0.0))
        for ak_a in aroma_keys:
            fluxes[f"v_aroma_{ak_a}"] = float(sol.fluxes.get(found_aromas[ak_a], 0.0))

    return fluxes, mode


# ══════════════════════════════════════════════════════════════════════════════════════
# HELPER 2 — RHS del sistema ODE
# ══════════════════════════════════════════════════════════════════════════════════════

def _compute_rhs_colloc(yk, fluxes, mode):
    """
    Evalua el RHS del sistema ODE dado el estado yk y los flujos FBA.
    Devuelve cdot (vector n_state, unidades fisicas).
    Ecuaciones identicas a celda 44 (_run_dfba_nrec).
    """
    cX    = max(0.0, yk[IDX_X])
    cProt = max(0.0, yk[IDX_PROT])
    cCarb = max(0.0, yk[IDX_CARB])
    cNrec = max(0.0, yk[IDX_NREC])
    Xa    = XA_FRACTION * cX

    rhs = np.zeros(n_state)

    rhs[IDX_X]     = max(0.0, fluxes["v_obj"]) * cX if mode == "BIOMASS" else 0.0
    rhs[IDX_NFREE] = -fluxes["vn_eff"] * cX if mode == "BIOMASS" else 0.0
    rhs[IDX_G]     = -MW_GLU * max(0.0, -fluxes["v_glu"]) * cX
    rhs[IDX_F]     = -MW_FRU * max(0.0, -fluxes["v_fru"]) * cX
    rhs[IDX_E]     =  MW_ETH * max(0.0,  fluxes["v_eth"]) * cX
    rhs[IDX_O2]    = -MW_O2  * max(0.0, -fluxes["v_o2"])  * cX
    rhs[IDX_PROT]  = (
        max(0.0, fluxes["v_obj"]) * cX * PROT_CONTENT_0
        - cProt * K_DEATH
        + Xa    * max(0.0, fluxes["v_prot"])
        - TURNOVER_LAMBDA * cProt
    )
    rhs[IDX_CARB]  = (
        max(0.0, fluxes["v_obj"]) * cX * CARB_CONTENT_0
        - cCarb * K_DEATH
    )

    nrec_inflow = Y_N_FROM_PROT * TURNOVER_LAMBDA * cProt
    if mode == "TURNOVER_ATPM":
        total_aa_mmol = sum(max(0.0, -fluxes.get(f"v_aa_{ak}", 0.0)) for ak in aa_keys)
        nrec_outflow  = total_aa_mmol * MW_N * cX
    else:
        nrec_outflow = 0.0
    rhs[IDX_NREC] = nrec_inflow - nrec_outflow

    for j_aa, ak in enumerate(aa_keys):
        v_aa = fluxes.get(f"v_aa_{ak}", 0.0)
        rhs[IDX_AA0 + j_aa] = cX * v_aa if mode == "BIOMASS" else 0.0

    for j_ar, ak_a in enumerate(aroma_keys):
        mw_a = AROMA_MW.get(ak_a, 0.1)
        rhs[IDX_AROMA0 + j_ar] = mw_a * max(0.0, fluxes.get(f"v_aroma_{ak_a}", 0.0)) * cX

    return rhs


# ══════════════════════════════════════════════════════════════════════════════════════
# FUNCION PRINCIPAL — dFBA SOA con colocaciones ortogonales Radau IIA
# ══════════════════════════════════════════════════════════════════════════════════════

def _run_dfba_colloc(
    *,
    scenario_name = "Anaerobico + N_rec + Radau IIA",
    aerobic_mode  = False,
    o2_init_g_l   = 0.0,
    nfe           = 18,    # elementos finitos (igual que NB2)
    ncp           = 3,     # puntos de colocacion Radau IIA (no cambiar)
    t_end         = 72.0,  # horizonte temporal [h]
    n_iter        = 2,     # iteraciones predictor-corrector por FE
                           #   1 = solo predictor (baja precision)
                           #   2 = 1 correccion (recomendado)
                           #   3+ = mayor precision, mas LPs
):
    """
    dFBA SOA con colocaciones ortogonales Radau IIA.

    ESTRUCTURA POR ELEMENTO FINITO i (i = 0..nfe-1)
    -------------------------------------------------
    c_prev : estado fisico al final del FE i-1 (o condicion inicial).
    h = t_end / nfe.
    Puntos temporales: t_loc[i,j] = t_start[i] + RADAU_ROOTS[j] * h.

    ITERACION PREDICTOR-CORRECTOR
    ------------------------------
      Prediccion:  C[:, i, j] = c_prev  para todo j.
      Para _iter en range(n_iter):
          Para j = 0..ncp-1:
              Resolver LP en C[:, i, j]  =>  V_flux[i][j], V_mode[i][j]
              CDOT[:, i, j] = RHS(C[:, i, j], V_flux[i][j])
          Actualizar (formula de colocacion Radau):
              C[:, i, :] = c_prev + h * CDOT[:, i, :] @ COLMAT_RADAU.T
      Endpoint del FE: c_prev <- C[:, i, ncp-1]

    CONSISTENCIA CDOT (diagnostico)
    ---------------------------------
      Por construccion:
          COLMAT_RADAU @ CDOT[:,i,:].T = (C[:,i,:] - c_prev) / h
      => cdot_colloc = COLMAT_INV @ (C - c_prev)/h identicamente igual a CDOT
      => ODE gap = 0 por construccion.
      Este warm-start no requiere saneamiento dinamico en NB2.

    SALIDAS
    --------
      states_df   : estados en todos los nfe*ncp+1 puntos de tiempo
      flux_df     : flujos en los nfe*ncp puntos de colocacion
      C           : (n_state, nfe, ncp) fisico  [warm_c para NB2]
      CDOT        : (n_state, nfe, ncp)          [warm_cdot para NB2]
      C_endpoints : (nfe+1, n_state) extremos de FEs
    """
    assert ncp == 3, "Esta implementacion asume ncp=3 (Radau IIA de 3 etapas)."

    h   = t_end / nfe
    hm  = np.full(nfe, h)
    nc  = n_state

    t_fe_start = np.concatenate([[0.0], np.cumsum(hm[:-1])])
    t_loc      = np.array([[t_fe_start[i] + RADAU_ROOTS[j] * h
                             for j in range(ncp)] for i in range(nfe)])

    C      = np.zeros((nc, nfe, ncp))
    CDOT   = np.zeros((nc, nfe, ncp))
    V_flux = [[None] * ncp for _ in range(nfe)]
    V_mode = [[""]   * ncp for _ in range(nfe)]

    N0_total   = float(globals().get("N0", 0.14))
    N0_ammonia = N0_total * N_AMMONIA_FRACTION
    N0_from_aa = max(0.0, N0_total - N0_ammonia)
    aa0_each   = (N0_from_aa / MW_N) / max(1, n_aa)

    # IC configurables via globales (hook para IC por experimento).
    # Defaults reproducen los valores historicos cuando los globales no estan seteados.
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

    c_prev      = c0.copy()
    C_endpoints = np.zeros((nfe + 1, nc))
    C_endpoints[0] = c0

    print(f"\n[{scenario_name}]  nfe={nfe}, ncp={ncp}, h={h:.3f} h, n_iter={n_iter}")
    print(f"  Total LPs (aprox) = {nfe} x {ncp} x {n_iter} = {nfe*ncp*n_iter}")

    for i in range(nfe):

        # Prediccion: todos los puntos de colocacion parten de c_prev
        for j in range(ncp):
            C[:, i, j] = np.maximum(c_prev, 0.0)

        for _it in range(n_iter):

            # 1. Resolver LP en cada punto de colocacion
            for j in range(ncp):
                # --- FASE 3: actualizar T_val segun T_PROFILE_DF o dynamic_temperature ---
                _t_cp = float(t_loc[i, j])
                _tp = globals().get("T_PROFILE_DF", None)
                if _tp is not None and len(_tp):
                    _T_now = float(np.interp(_t_cp,
                                             _tp["t_h"].to_numpy(),
                                             _tp["T_K"].to_numpy()))
                else:
                    _T_now = float(dynamic_temperature(_t_cp))
                globals()["T_val"] = _T_now

                fluxes, mode = _fba_at_state_colloc(C[:, i, j], aerobic_mode)

                if fluxes is None:
                    mode = "FAIL"
                    # Reusar flujos previos si existen; sino zero.
                    if V_flux[i][j] is not None:
                        fluxes = V_flux[i][j]
                    else:
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
                V_flux[i][j] = fluxes
                V_mode[i][j] = mode

                # 2. Derivada (RHS) en el punto de colocacion actual
                CDOT[:, i, j] = _compute_rhs_colloc(C[:, i, j], fluxes, mode)

            # 3. Actualizar estados via formula de colocacion Radau:
            #    C[s,i,j] = c_prev[s] + h * sum_k COLMAT[j,k]*CDOT[s,i,k]
            #    Vectorizado: (nc,ncp) = (nc,1) + h*(nc,ncp)@(ncp,ncp).T
            C_new = c_prev[:, None] + h * (CDOT[:, i, :] @ COLMAT_RADAU.T)
            C[:, i, :] = np.maximum(C_new, 0.0)

        # Endpoint del FE = ultimo punto de colocacion (tau=1 en Radau IIA)
        c_prev = C[:, i, ncp - 1].copy()

        # --- FASE 4: pulsos de N dentro de este FE ---
        _t_start_fe = h * i
        _t_end_fe   = h * (i + 1)
        _pulses = globals().get("PULSE_SCHEDULE", None) or []
        for _p in _pulses:
            _tp_pulse = float(_p["t_h"])
            if _t_start_fe < _tp_pulse <= _t_end_fe:
                _dn = float(_p.get("amount_gN_L", 0.0))
                # Suma al pool NFREE (ammonia). Distribucion entre AAs pendiente.
                c_prev[IDX_NFREE] += _dn

        C_endpoints[i + 1] = c_prev

    # Diagnostico ODE gap (debe ser 0 por construccion)
    max_ode_gap = 0.0
    for i in range(nfe):
        c_p = C_endpoints[i]
        for s in range(nc):
            rhs_j     = np.array([(C[s, i, j] - c_p[s]) / h for j in range(ncp)])
            cdot_coll = COLMAT_RADAU_INV @ rhs_j
            gap = float(np.max(np.abs(cdot_coll - CDOT[s, i, :])))
            if gap > max_ode_gap:
                max_ode_gap = gap

    # ══════════════════════════════════════════════════════════════════════════════
    # DESEMPAQUETAR RESULTADOS
    # ══════════════════════════════════════════════════════════════════════════════

    t_all  = np.concatenate([[0.0], t_loc.ravel()])
    n_pts  = 1 + nfe * ncp
    Y_all  = np.zeros((n_pts, nc))
    Y_all[0] = c0
    for i in range(nfe):
        for j in range(ncp):
            Y_all[1 + i * ncp + j] = C[:, i, j]

    t_colloc = t_loc.ravel()

    def _flux_series(key):
        out = np.zeros(nfe * ncp)
        for i in range(nfe):
            for j in range(ncp):
                f = V_flux[i][j]
                out[i * ncp + j] = f.get(key, 0.0) if f else 0.0
        return out

    mu_c    = _flux_series("v_obj")
    vglu_c  = np.maximum(0.0, -_flux_series("v_glu"))
    vfru_c  = np.maximum(0.0, -_flux_series("v_fru"))
    veth_c  = np.maximum(0.0,  _flux_series("v_eth"))
    vo2_c   = np.maximum(0.0, -_flux_series("v_o2"))
    vn_c    = np.maximum(0.0,  _flux_series("vn_eff"))
    vatpm_c = np.maximum(0.0,  _flux_series("v_atpm"))
    vprot_c = np.maximum(0.0,  _flux_series("v_prot"))
    gam_c   = _flux_series("full_gam")
    mode_c  = [V_mode[i][j] for i in range(nfe) for j in range(ncp)]

    aa_uptake_c  = {ak: np.maximum(0.0, -_flux_series(f"v_aa_{ak}"))    for ak in aa_keys}
    aroma_flux_c = {ak: np.maximum(0.0,  _flux_series(f"v_aroma_{ak}")) for ak in aroma_keys}

    stat_mask = np.array([m == "TURNOVER_ATPM" for m in mode_c])
    t_stat    = float(t_colloc[np.argmax(stat_mask)]) if stat_mask.any() else float("nan")

    X_all     = Y_all[:, IDX_X]
    Nfree_all = Y_all[:, IDX_NFREE]
    Prot_all  = Y_all[:, IDX_PROT]
    Carb_all  = Y_all[:, IDX_CARB]
    Nrec_all  = Y_all[:, IDX_NREC]

    aa_conc_all  = {ak: Y_all[:, IDX_AA0 + j]     for j,   ak  in enumerate(aa_keys)}
    aroma_all    = {ak: Y_all[:, IDX_AROMA0 + j_a] for j_a, ak in enumerate(aroma_keys)}

    aa_n_all        = (MW_N * np.sum(np.column_stack([aa_conc_all[ak] for ak in aa_keys]), axis=1)
                       if aa_keys else np.zeros(n_pts))
    N_in_prot_all   = Y_N_FROM_PROT * Prot_all
    N_total_sys_all = Nfree_all + aa_n_all + N_in_prot_all + Nrec_all
    N_bal_err       = abs(float(N_total_sys_all[-1]) - float(N_total_sys_all[0]))

    P_frac_all = Prot_all / np.maximum(X_all, EPS_C)
    C_frac_all = Carb_all / np.maximum(X_all, EPS_C)
    gam_full_all = np.array([
        _compute_full_gam(float(P_frac_all[ii]), RNA_FRAC, float(C_frac_all[ii]),
                          Pbase_global, Rbase_global, Cbase_global)
        for ii in range(n_pts)
    ])

    biomass_steps  = sum(m == "BIOMASS"       for m in mode_c)
    turnover_steps = sum(m == "TURNOVER_ATPM" for m in mode_c)

    print(f"  Puntos col: BIOMASS={biomass_steps}, TURNOVER_ATPM={turnover_steps}")
    print(f"  t_stationary ~ {t_stat:.1f} h")
    print(f"  X_final={X_all[-1]:.4f} gDW/L  |  E_final={Y_all[-1, IDX_E]:.4f} g/L")
    print(f"  Prot_final={Prot_all[-1]:.4f} g/L  |  Prot/X={P_frac_all[-1]:.4f}")
    print(f"  N_rec_final={Nrec_all[-1]:.6f} gN/L")
    print(f"  GAM: inicial={gam_full_all[0]:.3f} -> final={gam_full_all[-1]:.3f} mmol ATP/gDW")
    print(f"  Balance N: {N_total_sys_all[0]:.6f} -> {N_total_sys_all[-1]:.6f} gN/L  "
          f"(err={N_bal_err:.2e}  {'OK' if N_bal_err < 0.01 else 'CHECK!'})")
    print(f"  ODE gap max = {max_ode_gap:.3e}  "
          f"({'OK (construccion)' if max_ode_gap < 1e-10 else 'CHECK'})")

    # ══════════════════════════════════════════════════════════════════════════════
    # DataFrames
    # ══════════════════════════════════════════════════════════════════════════════

    states_df = pd.DataFrame({
        "t_h":                 t_all,
        "X_gDW_L":             X_all,
        "N_gN_L":              Nfree_all + aa_n_all,
        "N_free_gN_L":         Nfree_all,
        "N_from_AA_gN_L":      aa_n_all,
        "N_rec_gN_L":          Nrec_all,
        "N_in_Prot_gN_L":      N_in_prot_all,
        "N_total_system_gN_L": N_total_sys_all,
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
        "t_h":                 t_colloc,
        "fe_idx":              np.repeat(np.arange(nfe), ncp),
        "cp_idx":              np.tile(np.arange(ncp), nfe),
        "mode":                mode_c,
        "mu_h_inv":            mu_c,
        "v_glu":               vglu_c,
        "v_fru":               vfru_c,
        "v_eth":               veth_c,
        "v_o2":                vo2_c,
        "v_N_uptake_gN_gDW_h": vn_c,
        "v_ATPM":              vatpm_c,
        "v_Prot":              vprot_c,
        "GAM":                 gam_c,
    })
    for ak in aa_keys:
        flux_df[f"v_uptake_{ak}"] = aa_uptake_c[ak]
    for ak in aroma_keys:
        flux_df[f"v_{ak}"] = aroma_flux_c[ak]

    # Fluxes extra por rid (acetate esters / ethyl acetate): exporta v_aroma_<rid>.
    for _rid in globals().get("EXTRA_AROMA_RIDS", []):
        flux_df[f"v_aroma_{_rid}"] = np.maximum(0.0, _flux_series(f"v_aroma_{_rid}"))

    display(flux_df.head(9))
    display(flux_df.tail(9))

    # ══════════════════════════════════════════════════════════════════════════════
    # GRAFICOS
    # ══════════════════════════════════════════════════════════════════════════════

    state_series = [
        ("X (gDW/L)",      X_all,              np.interp(t_all, t_colloc, mu_c),    "mu (1/h)"),
        ("N_ext (gN/L)",   Nfree_all+aa_n_all, np.interp(t_all, t_colloc, vn_c),    "v_N"),
        ("N_rec (gN/L)",   Nrec_all,           np.zeros(n_pts),                      "None"),
        ("Glucosa (g/L)",  Y_all[:,IDX_G],     np.interp(t_all, t_colloc, vglu_c),  "v_glu"),
        ("Fructosa (g/L)", Y_all[:,IDX_F],     np.interp(t_all, t_colloc, vfru_c),  "v_fru"),
        ("Etanol (g/L)",   Y_all[:,IDX_E],     np.interp(t_all, t_colloc, veth_c),  "v_eth"),
        ("O2 (g/L)",       Y_all[:,IDX_O2],    np.interp(t_all, t_colloc, vo2_c),   "v_o2"),
        ("Prot (g/L)",     Prot_all,           np.interp(t_all, t_colloc, vprot_c), "v_Prot"),
        ("Carb (g/L)",     Carb_all,           np.zeros(n_pts),                      "None"),
        ("N total (gN/L)", N_total_sys_all,    np.zeros(n_pts),                      "balance"),
    ] + [
        (f"AA_{ak} (mmol/L)", aa_conc_all[ak], aa_uptake_c[ak], f"v_{ak}")
        for ak in aa_keys
    ] + [
        (f"{ak} (mg/L)", 1000.0*aroma_all[ak], aroma_flux_c[ak], f"v_{ak}")
        for ak in aroma_keys
    ]

    ncols_s = 4
    nrows_s = math.ceil(len(state_series) / ncols_s)
    fig, ax  = plt.subplots(nrows_s, ncols_s,
                             figsize=(5.2*ncols_s, 3.5*nrows_s), sharex=True, dpi=120)
    ax = np.atleast_1d(ax).ravel()
    for i_p, (ttl, yy_, vv_, vlbl) in enumerate(state_series):
        ax[i_p].plot(t_all, yy_, lw=2, color="tab:blue")
        ax[i_p].scatter(t_colloc, np.interp(t_colloc, t_all, yy_),
                         s=14, color="tab:green", zorder=5, alpha=0.75)
        if np.isfinite(t_stat):
            ax[i_p].axvline(t_stat, ls="--", lw=1, color="navy", alpha=0.5)
        ax[i_p].set_title(ttl, fontsize=9)
        ax[i_p].grid(alpha=0.3)
        ax[i_p].set_xlabel("t (h)", fontsize=8)
        ax2 = ax[i_p].twinx()
        vv_cp = vv_ if len(vv_) == len(t_colloc) else np.interp(t_colloc, t_all, vv_)
        ax2.step(t_colloc, vv_cp, where="post", lw=1.4, ls="--",
                 color="tab:orange", alpha=0.85)
        ax2.set_ylabel(vlbl, color="tab:orange", fontsize=7)
        ax2.tick_params(axis="y", labelcolor="tab:orange", labelsize=7)
    for j_p in range(i_p + 1, len(ax)):
        ax[j_p].axis("off")
    plt.suptitle(
        f"dFBA Radau IIA [{nfe} FE, {ncp} CP, {n_iter} iter] — {scenario_name}",
        fontsize=11, y=1.01)
    plt.tight_layout()
    plt.show()

    # GAM + composicion
    fig_g, ax_g = plt.subplots(1, 3, figsize=(15, 4), dpi=120)
    ax_g[0].plot(t_all, gam_full_all, lw=2.5, color="tab:red")
    ax_g[0].axhline(30.49, ls="--", lw=1.2, color="gray", label="GAM_base")
    if np.isfinite(t_stat):
        ax_g[0].axvline(t_stat, ls=":", lw=1.5, color="navy", alpha=0.6)
    ax_g[0].set_title("GAM (mmol ATP/gDW)")
    ax_g[0].set_xlabel("t (h)")
    ax_g[0].legend(fontsize=8)
    ax_g[0].grid(alpha=0.3)

    ax_g[1].plot(t_all, P_frac_all, lw=2.5, color="tab:green", label="Protein (g/gDW)")
    ax_g[1].plot(t_all, C_frac_all, lw=2.5, color="tab:purple", label="Carb (g/gDW)")
    ax_g[1].axhline(Pbase_global, ls="--", lw=1, color="tab:green", alpha=0.5)
    ax_g[1].axhline(Cbase_global, ls="--", lw=1, color="tab:purple", alpha=0.5)
    if np.isfinite(t_stat):
        ax_g[1].axvline(t_stat, ls=":", lw=1.5, color="navy", alpha=0.6)
    ax_g[1].set_title("Composicion celular (g/gDW)")
    ax_g[1].set_xlabel("t (h)")
    ax_g[1].legend(fontsize=8)
    ax_g[1].grid(alpha=0.3)

    ax_g[2].plot(t_all, Prot_all, lw=2.5, color="tab:green", label="Prot (g/L)")
    ax_g[2].plot(t_all, X_all,    lw=2.5, color="tab:blue",  label="X (gDW/L)")
    ax_g[2].plot(t_all, Nrec_all, lw=2.5, color="tab:red",   label="N_rec (gN/L)")
    if np.isfinite(t_stat):
        ax_g[2].axvline(t_stat, ls=":", lw=1.5, color="navy", alpha=0.6)
    ax_g[2].set_title("Biomasa, Proteina y N_rec")
    ax_g[2].set_xlabel("t (h)")
    ax_g[2].legend(fontsize=8)
    ax_g[2].grid(alpha=0.3)
    plt.suptitle(f"GAM variable – {scenario_name}", fontsize=11)
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
            if np.isfinite(t_stat):
                ax_a[i_a].axvline(t_stat, ls=":", lw=1.2, color="navy", alpha=0.5)
        for j_a in range(i_a + 1, len(ax_a)):
            ax_a[j_a].axis("off")
        plt.suptitle(f"Aromas – {scenario_name}", fontsize=11)
        plt.tight_layout()
        plt.show()

    # Guardar CSV
    sfx   = "aerobic_colloc" if aerobic_mode else "anaerobic_colloc"
    out_s = OUT_DIR / f"dfba_colloc_states_{sfx}.csv"
    out_f = OUT_DIR / f"dfba_colloc_flux_{sfx}.csv"
    states_df.to_csv(out_s, index=False)
    flux_df.to_csv(out_f,   index=False)
    print(f"\nArchivos guardados:\n  {out_s}\n  {out_f}")

    return dict(
        states_df   = states_df,
        flux_df     = flux_df,
        C           = C,
        CDOT        = CDOT,
        C_endpoints = C_endpoints,
        V_flux      = V_flux,
        V_mode      = V_mode,
        t_loc       = t_loc,
        hm          = hm,
        t_stat      = t_stat,
        nfe         = nfe,
        ncp         = ncp,
        max_ode_gap = max_ode_gap,
    )


# ══════════════════════════════════════════════════════════════════════════════════════
# EJECUCION
# ══════════════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("Escenario: Anaerobico + N_rec + Radau IIA  (nfe=18, ncp=3, n_iter=2)")
print("=" * 80)

res_colloc_ana = _run_dfba_colloc(
    scenario_name = "Anaerobico + N_rec + Radau IIA",
    aerobic_mode  = False,
    o2_init_g_l   = 0.0,
    nfe           = 18,
    ncp           = 3,
    t_end         = 72.0,
    n_iter        = 2,
)

# ══════════════════════════════════════════════════════════════════════════════════════
# EXPORTAR WARM-START EN FORMATO NB2
# ══════════════════════════════════════════════════════════════════════════════════════
# C[n_state, nfe, ncp] y CDOT[n_state, nfe, ncp] satisfacen por construccion:
#   COLMAT_RADAU @ CDOT[:,i,:].T = (C[:,i,:] - C_endpoints[i,:]) / h
# => ODE gap = 0 => no se requiere saneamiento dinamico adicional en NB2.
# Formato CSV aplanado (columnas: state_idx, fe, cp, value).

_C_ws    = res_colloc_ana["C"]
_CDOT_ws = res_colloc_ana["CDOT"]
_hm_ws   = res_colloc_ana["hm"]
_nfe_ws  = res_colloc_ana["nfe"]
_ncp_ws  = res_colloc_ana["ncp"]

_rows_c    = []
_rows_cdot = []
for _s in range(n_state):
    for _i in range(_nfe_ws):
        for _j in range(_ncp_ws):
            _rows_c.append(   {"state_idx": _s, "fe": _i, "cp": _j, "value": float(_C_ws[_s, _i, _j])})
            _rows_cdot.append({"state_idx": _s, "fe": _i, "cp": _j, "value": float(_CDOT_ws[_s, _i, _j])})

_ws_c_file    = OUT_DIR / "warm_start_states_robust_colloc.csv"
_ws_cdot_file = OUT_DIR / "warm_start_cdot_robust_colloc.csv"
_ws_h_file    = OUT_DIR / "warm_start_h_robust_fe.csv"

pd.DataFrame(_rows_c).to_csv(_ws_c_file,    index=False)
pd.DataFrame(_rows_cdot).to_csv(_ws_cdot_file, index=False)
np.savetxt(_ws_h_file, _hm_ws, delimiter=",")

print(f"\nWarm-start Radau IIA exportado (listo para NB2 sin saneamiento dinamico):")
print(f"  {_ws_c_file}")
print(f"  {_ws_cdot_file}")
print(f"  {_ws_h_file}")
print(f"  ODE gap = {res_colloc_ana['max_ode_gap']:.2e}")
