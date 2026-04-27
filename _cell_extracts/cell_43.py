import math
import re
import numpy as np
import pandas as pd
import cobra

# ======================================================================================
# dFBA SOA + pFBA con protein-turnover CONSERVATIVO (pool N_rec) + aromas
# ======================================================================================

import matplotlib.pyplot as plt

print("\n" + "=" * 100)
print("dFBA SOA + pFBA con GAM variable + protein-turnover CONSERVATIVO (N_rec) + aromas")
print("=" * 100)

# ── Verificar dependencias de celdas anteriores ────────────────────────────────────
_required_sym_nrec = [
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
]
_missing_nrec = [s for s in _required_sym_nrec if s not in globals()]
if _missing_nrec:
    raise RuntimeError(f"Faltan símbolos requeridos: {_missing_nrec}")

_kinetic_limits_fn = globals().get("_kinetic_limits_fn", globals().get("_kinetic_limits"))
N_VITAMIN_EX_IDS_LOC = globals().get("N_VITAMIN_EX_IDS_LOC", globals().get("N_VITAMIN_EX_IDS", []))
EPS_LOC = float(globals().get("EPS_LOC", globals().get("EPS", 1e-9)))

# ══════════════════════════════════════════════════════════════════════════════════════
# PARÁMETROS
# ══════════════════════════════════════════════════════════════════════════════════════

# --- Composición y turnover ---
PROT_CONTENT_0   = 0.46    # g protein / gDW  (initial)
CARB_CONTENT_0   = 0.37    # g carbohydrate / gDW  (initial)
K_DEATH          = 0.005   # 1/h  protein decay
TURNOVER_LAMBDA  = 0.03    # 1/h  protein turnover rate
VPROT_FLOOR_FRAC = 0.90    # fraction of max v_Prot to lock as floor
XA_FRACTION      = 1.0     # active biomass fraction

# --- Nitrogen partitioning ---
N_AMMONIA_FRACTION          = 0.50   # fraction of total N0 as free ammonia
N_TOTAL_DEPLETION_THRESHOLD = 1e-3   # gN/L  below this -> stationary phase

# --- AA uptake during growth ---
K_AA_UPTAKE_GROWTH = 0.08   # 1/h  kinetic uptake of AA during growth
BIOMASS_LOCK_FRAC  = 0.999  # lock biomass flux for phase B
AA_UPTAKE_WEIGHT   = 1.0    # weight for secondary AA-uptake objective

# --- NEW: Recycled nitrogen pool (N_rec) ---
Y_N_FROM_PROT   = 0.16     # gN / gProtein  (avg amino acid N content)
K_NREC_UPTAKE   = 0.10     # 1/h  max specific rate of N_rec consumption
N_REC_INIT      = 0.0      # gN/L  initial recycled nitrogen

# --- AA alpha (fractional allocation per precursor) ---
aa_keys   = list(found_precursors.keys())
aroma_keys = list(found_aromas.keys())
n_aa      = len(aa_keys)
n_aroma   = len(aroma_keys)

AA_ALPHA_STATIC = 1.0 / max(1, n_aa)
AA_ALPHA        = {k: AA_ALPHA_STATIC for k in aa_keys}

# --- Build tracked reaction IDs ---
_tracked_ids = list(dict.fromkeys(
    [OBJ_ID, GLU_ID, FRU_ID, ETH_ID, O2_ID, ATPM_ID, PROT_RXN_ID]
    + list(found_precursors.values())
    + list(found_aromas.values())
    + list(KINETIC_N_SOURCE_IDS)
))
_tracked_ids = [rid for rid in _tracked_ids if rid in rxn_by_id]

_base_lb_nrec = {rid: float(rxn_by_id[rid].lower_bound) for rid in _tracked_ids}
_base_ub_nrec = {rid: float(rxn_by_id[rid].upper_bound) for rid in _tracked_ids}

print(f"TURNOVER_LAMBDA      = {TURNOVER_LAMBDA}")
print(f"Y_N_FROM_PROT        = {Y_N_FROM_PROT} gN/gProt")
print(f"K_NREC_UPTAKE        = {K_NREC_UPTAKE} 1/h")
print(f"N_REC_INIT           = {N_REC_INIT} gN/L")
print(f"N_AMMONIA_FRACTION   = {N_AMMONIA_FRACTION}")
print(f"Precursors tracked   = {n_aa}  |  Aromas tracked = {n_aroma}")

# ══════════════════════════════════════════════════════════════════════════════════════
# MAIN dFBA FUNCTION
# ══════════════════════════════════════════════════════════════════════════════════════

def _run_dfba_nrec(
    *,
    scenario_name="Anaeróbico + N_rec conservativo",
    aerobic_mode=False,
    o2_init_g_l=0.0,
    dfba_dt=1.0,
    t_end=72.0,
):
    """
    dFBA with conservative protein turnover via explicit N_rec pool.

    State vector layout:
      [X, N_free, G, F, E, O2, Prot, Carb, N_rec, AA_0..AA_{n_aa-1}, Aroma_0..Aroma_{n_aroma-1}]
       0    1     2  3  4   5    6     7      8     9..               9+n_aa..

    Key change vs. previous version:
      - Protein degradation (λ·Prot) feeds N_rec, not directly AA uptake caps.
      - In stationary phase, AA uptake caps are set by N_rec availability:
            q_total_nrec = K_NREC_UPTAKE * N_rec / X
        distributed uniformly among precursors (or via AA_ALPHA).
      - N_rec is consumed when the FBA actually uses those AA fluxes.
      - Mass balance: N leaving Prot = N entering N_rec = N consumed by metabolism.
    """
    dt   = float(dfba_dt)
    tt   = np.arange(0.0, t_end + 0.5 * dt, dt)
    n_st = len(tt) - 1

    N0_total   = float(globals().get("N0", 0.14))
    N0_ammonia = N0_total * N_AMMONIA_FRACTION
    N0_from_aa = max(0.0, N0_total - N0_ammonia)
    aa0_each   = (N0_from_aa / MW_N) / max(1, n_aa)

    # State vector dimensions
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

    Y = np.zeros((len(tt), n_state), dtype=float)
    Y[0, IDX_X]     = 0.5
    Y[0, IDX_NFREE] = N0_ammonia
    Y[0, IDX_G]     = 110.0
    Y[0, IDX_F]     = 110.0
    Y[0, IDX_E]     = 0.0
    Y[0, IDX_O2]    = o2_init_g_l
    Y[0, IDX_PROT]  = Y[0, IDX_X] * PROT_CONTENT_0
    Y[0, IDX_CARB]  = Y[0, IDX_X] * CARB_CONTENT_0
    Y[0, IDX_NREC]  = N_REC_INIT
    for j in range(n_aa):
        Y[0, IDX_AA0 + j] = aa0_each

    # History arrays
    mu_h     = np.zeros(n_st)
    vglu_h   = np.zeros(n_st)
    vfru_h   = np.zeros(n_st)
    veth_h   = np.zeros(n_st)
    vo2_h    = np.zeros(n_st)
    vn_h     = np.zeros(n_st)
    vatpm_h  = np.zeros(n_st)
    vprot_h  = np.zeros(n_st)
    gam_h    = np.zeros(n_st)
    p_frac_h = np.zeros(n_st)
    c_frac_h = np.zeros(n_st)

    # N_rec-specific histories
    nrec_inflow_h   = np.zeros(n_st)   # gN/L/h entering N_rec from turnover
    nrec_outflow_h  = np.zeros(n_st)   # gN/L/h leaving N_rec to metabolism
    nrec_cap_h      = np.zeros(n_st)   # max allowed v_Nrec (mmol equiv / gDW / h)

    aroma_fh    = {k: np.zeros(n_st) for k in aroma_keys}
    aa_uptake_h = {k: np.zeros(n_st) for k in aa_keys}
    aa_net_h    = {k: np.zeros(n_st) for k in aa_keys}
    mode_h      = []

    for k in range(n_st):
        yk = np.maximum(Y[k, :], 0.0)
        cX     = yk[IDX_X]
        cN_free = yk[IDX_NFREE]
        cG     = yk[IDX_G]
        cF     = yk[IDX_F]
        cE     = yk[IDX_E]
        cO2    = yk[IDX_O2]
        cProt  = yk[IDX_PROT]
        cCarb  = yk[IDX_CARB]
        cNrec  = yk[IDX_NREC]
        aa_state = {ak: float(yk[IDX_AA0 + j]) for j, ak in enumerate(aa_keys)}

        # Total external N (free ammonia + AA in medium)
        aa_n_sum = MW_N * sum(aa_state.values())
        cN_total = cN_free + aa_n_sum

        # Dynamic protein/carb fractions
        P_frac = cProt / max(cX, EPS_LOC)
        C_frac = cCarb / max(cX, EPS_LOC)
        R_frac = RNA_FRAC

        full_gam = _compute_full_gam(P_frac, R_frac, C_frac,
                                      Pbase_global, Rbase_global, Cbase_global)
        gam_h[k]    = full_gam
        p_frac_h[k] = P_frac
        c_frac_h[k] = C_frac

        lim_glu, lim_fru, lim_eth, lim_obj, lim_n_map = _kinetic_limits_fn(
            cX, cN_free, cG, cF, cE, cO2)

        with model as mtmp:
            rb = {r.id: r for r in mtmp.reactions}

            # Reset tracked bounds
            for rid in _tracked_ids:
                if rid in rb:
                    rb[rid].lower_bound = _base_lb_nrec[rid]
                    rb[rid].upper_bound = _base_ub_nrec[rid]

            # Sugar kinetics
            if GLU_ID in rb:
                rb[GLU_ID].lower_bound = max(rb[GLU_ID].lower_bound, -lim_glu)
                rb[GLU_ID].upper_bound = min(rb[GLU_ID].upper_bound, 0.0)
            if FRU_ID in rb:
                rb[FRU_ID].lower_bound = max(rb[FRU_ID].lower_bound, -lim_fru)
                rb[FRU_ID].upper_bound = min(rb[FRU_ID].upper_bound, 0.0)

            # N kinetic limits on exchange reactions
            for rid in KINETIC_N_SOURCE_IDS:
                if rid in rb:
                    Li = float(lim_n_map.get(rid, 0.0))
                    rb[rid].lower_bound = max(rb[rid].lower_bound, -Li)
                    rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)

            # O2
            if aerobic_mode:
                lim_o2 = O2_VMAX_UPTAKE * (cO2 / (cO2 + O2_KO_G_L + EPS_LOC))
                lim_o2 = 0.0 if cO2 <= O2_DEPLETION_THRESHOLD else max(0.0, lim_o2)
                rb[O2_ID].lower_bound = -lim_o2
                rb[O2_ID].upper_bound = 0.0
            else:
                rb[O2_ID].lower_bound = 0.0
                rb[O2_ID].upper_bound = 0.0

            # Product caps
            if APPLY_PRODUCT_CAPS:
                rb[ETH_ID].upper_bound = min(rb[ETH_ID].upper_bound, lim_eth)
                rb[OBJ_ID].upper_bound = min(rb[OBJ_ID].upper_bound, lim_obj)
            else:
                rb[ETH_ID].upper_bound = max(rb[ETH_ID].upper_bound, 1000.0)
                rb[OBJ_ID].upper_bound = max(rb[OBJ_ID].upper_bound, 1000.0)

            # Apply variable GAM
            _apply_variable_gam(mtmp, full_gam)

            # ══════════════════════════════════════════════════════════
            # PHASE SWITCH + AA UPTAKE CAPS
            # ══════════════════════════════════════════════════════════

            if cN_total > N_TOTAL_DEPLETION_THRESHOLD:
                # ── GROWTH PHASE (BIOMASS) ──
                # AA uptake from extracellular pool: kinetic Michaelis-like
                for ak, rid in found_precursors.items():
                    if rid in rb:
                        q_i = K_AA_UPTAKE_GROWTH * aa_state[ak] / max(cX, EPS_LOC)
                        q_i = max(0.0, q_i)
                        rb[rid].lower_bound = -q_i
                        rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)

                mode_h.append("BIOMASS")
                mtmp.objective = OBJ_ID
                mtmp.objective_direction = "max"
                try:
                    sol_mu = mtmp.optimize()
                except Exception:
                    mode_h[-1] = "FAIL"; continue
                if getattr(sol_mu, "status", "") != "optimal":
                    mode_h[-1] = "NO_OPT"; continue

                mu_star = max(0.0, float(sol_mu.fluxes.get(OBJ_ID, 0.0)))
                rb[OBJ_ID].lower_bound = max(rb[OBJ_ID].lower_bound,
                                              BIOMASS_LOCK_FRAC * mu_star)

                # Secondary: maximize AA uptake at locked growth
                uptake_obj = {mtmp.reactions.get_by_id(rid): -AA_UPTAKE_WEIGHT
                              for rid in found_precursors.values() if rid in rb}
                mtmp.objective = uptake_obj
                mtmp.objective_direction = "max"
                try:
                    sol = cobra.flux_analysis.pfba(mtmp, fraction_of_optimum=1.0)
                except Exception:
                    mode_h[-1] = "FAIL"; continue
                if getattr(sol, "status", "optimal") != "optimal":
                    mode_h[-1] = "NO_OPT"; continue

                nrec_cap_h[k] = 0.0  # not used in growth

            else:
                # ── STATIONARY PHASE (TURNOVER_ATPM) ──
                # Close all external N sources
                for rid in KINETIC_N_SOURCE_IDS:
                    if rid in rb:
                        rb[rid].lower_bound = max(rb[rid].lower_bound, 0.0)
                        rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)
                for rid in N_VITAMIN_EX_IDS_LOC:
                    if rid in rb:
                        rb[rid].lower_bound = max(rb[rid].lower_bound, 0.0)
                        rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)

                # No growth
                rb[OBJ_ID].lower_bound = 0.0
                rb[OBJ_ID].upper_bound = 0.0

                # ── KEY CHANGE: AA uptake caps from N_rec pool ──
                # Total specific N_rec available for AA uptake:
                #   q_Nrec_total = K_NREC_UPTAKE * N_rec / X   [gN / gDW / h]
                # Convert to mmol AA / gDW / h assuming average MW_N per AA:
                #   q_AA_total = q_Nrec_total / MW_N  [mmol_N / gDW / h]
                # Distribute among precursors via AA_ALPHA
                q_Nrec_total = K_NREC_UPTAKE * cNrec / max(cX, EPS_LOC)  # gN/gDW/h
                q_Nrec_total = max(0.0, q_Nrec_total)
                q_AA_total_mmol = q_Nrec_total / MW_N  # mmol_N/gDW/h ≈ mmol_AA/gDW/h

                nrec_cap_h[k] = q_Nrec_total

                for ak, rid in found_precursors.items():
                    if rid in rb:
                        q_i = q_AA_total_mmol * AA_ALPHA[ak]
                        q_i = max(0.0, q_i)
                        rb[rid].lower_bound = -q_i
                        rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)

                # Maximize protein resynthesis first
                vprot_star = 0.0
                try:
                    mtmp.objective = PROT_RXN_ID
                    mtmp.objective_direction = "max"
                    sol_prot = mtmp.optimize()
                    if getattr(sol_prot, "status", "") == "optimal":
                        vprot_star = max(0.0, float(sol_prot.fluxes.get(PROT_RXN_ID, 0.0)))
                except Exception:
                    vprot_star = 0.0

                if vprot_star > 1e-10:
                    rb[PROT_RXN_ID].lower_bound = max(
                        rb[PROT_RXN_ID].lower_bound,
                        VPROT_FLOOR_FRAC * vprot_star
                    )

                # ATPM bounds for stationary
                rb[ATPM_ID].lower_bound = max(rb[ATPM_ID].lower_bound, ATPM_LB_NO_GROWTH)
                rb[ATPM_ID].upper_bound = ATPM_UB_NO_GROWTH

                mtmp.objective = ATPM_ID
                mtmp.objective_direction = "max"
                mode_h.append("TURNOVER_ATPM")

                try:
                    sol = cobra.flux_analysis.pfba(mtmp, fraction_of_optimum=1.0)
                except Exception:
                    mode_h[-1] = "FAIL"; continue
                if getattr(sol, "status", "optimal") != "optimal":
                    mode_h[-1] = "NO_OPT"; continue

            # ══════════════════════════════════════════════════════════
            # EXTRACT FLUXES
            # ══════════════════════════════════════════════════════════

            v_obj  = float(sol.fluxes.get(OBJ_ID, 0.0))
            v_glu  = float(sol.fluxes.get(GLU_ID, 0.0))
            v_fru  = float(sol.fluxes.get(FRU_ID, 0.0))
            v_eth  = float(sol.fluxes.get(ETH_ID, 0.0))
            v_o2   = float(sol.fluxes.get(O2_ID, 0.0))
            v_atpm = float(sol.fluxes.get(ATPM_ID, 0.0))
            v_prot = float(sol.fluxes.get(PROT_RXN_ID, 0.0))

            # N uptake from external sources (growth phase only)
            if mode_h[-1] == "BIOMASS":
                vn_eff = 0.0
                for rid in KINETIC_N_SOURCE_IDS:
                    v_i = float(sol.fluxes.get(rid, 0.0))
                    vn_eff += max(0.0, -v_i) * N_atoms_map.get(rid, 1.0) * MW_N
            else:
                vn_eff = 0.0

            # AA precursor fluxes
            total_aa_uptake_mmol = 0.0
            for j_aa, ak in enumerate(aa_keys):
                rid = found_precursors[ak]
                v_aa = float(sol.fluxes.get(rid, 0.0))
                aa_uptake_h[ak][k] = max(0.0, -v_aa)
                aa_net_h[ak][k]    = v_aa
                total_aa_uptake_mmol += max(0.0, -v_aa)

            # Aroma fluxes
            for ak_a in aroma_keys:
                rid_a = found_aromas[ak_a]
                aroma_fh[ak_a][k] = max(0.0, float(sol.fluxes.get(rid_a, 0.0)))

        # ══════════════════════════════════════════════════════════════
        # STORE FLUX HISTORIES
        # ══════════════════════════════════════════════════════════════

        mu_h[k]    = max(0.0, v_obj)
        vglu_h[k]  = max(0.0, -v_glu)
        vfru_h[k]  = max(0.0, -v_fru)
        veth_h[k]  = max(0.0, v_eth)
        vo2_h[k]   = max(0.0, -v_o2)
        vn_h[k]    = max(0.0, vn_eff)
        vatpm_h[k] = max(0.0, v_atpm)
        vprot_h[k] = max(0.0, v_prot)

        # ══════════════════════════════════════════════════════════════
        # ODE INTEGRATION (Euler)
        # ══════════════════════════════════════════════════════════════

        # --- Biomass ---
        dX = v_obj * cX if mode_h[-1] == "BIOMASS" else 0.0

        # --- Free nitrogen ---
        dN_free = -vn_eff * cX if mode_h[-1] == "BIOMASS" else 0.0

        # --- Sugars/Ethanol/O2 ---
        dG  = -MW_GLU * max(0.0, -v_glu) * cX
        dF  = -MW_FRU * max(0.0, -v_fru) * cX
        dE  =  MW_ETH * max(0.0, v_eth) * cX
        dO2 = -MW_O2  * max(0.0, -v_o2) * cX

        # --- Protein (unchanged equation) ---
        Xa = XA_FRACTION * cX
        dProt = (
            max(0.0, v_obj) * cX * PROT_CONTENT_0   # new protein from growth
            - cProt * K_DEATH                         # death/decay
            + Xa * max(0.0, v_prot)                   # resynthesis from recycled AA
            - TURNOVER_LAMBDA * cProt                  # degradation -> feeds N_rec
        )

        # --- Carbohydrate ---
        dCarb = max(0.0, v_obj) * cX * CARB_CONTENT_0 - cCarb * K_DEATH

        # --- N_rec: recycled nitrogen pool ---
        # Inflow: degraded protein releases nitrogen
        nrec_inflow = Y_N_FROM_PROT * TURNOVER_LAMBDA * cProt  # gN/L/h

        # Outflow: metabolism consumes recycled N via AA uptake
        # In stationary: actual AA uptake by FBA * MW_N converts to gN consumed
        # In growth: N_rec is not the source, so outflow = 0
        if mode_h[-1] == "TURNOVER_ATPM":
            # total_aa_uptake_mmol is in mmol/gDW/h; multiply by X for volumetric
            # Each mmol AA ≈ 1 mmol N (simplification; ~1 N atom per AA on avg)
            nrec_outflow = total_aa_uptake_mmol * MW_N * cX  # gN/L/h
        else:
            nrec_outflow = 0.0

        nrec_inflow_h[k]  = nrec_inflow
        nrec_outflow_h[k] = nrec_outflow

        dNrec = nrec_inflow - nrec_outflow

        # --- Extracellular amino acids ---
        # Growth phase: AA consumed from medium
        # Stationary: AA exchanges are closed for external; only N_rec feeds FBA
        #   -> no change in extracellular AA pools during stationary

        # --- Integration ---
        Y[k+1, IDX_X]     = max(0.0, cX     + dt * dX)
        Y[k+1, IDX_NFREE] = max(0.0, cN_free + dt * dN_free)
        Y[k+1, IDX_G]     = max(0.0, cG     + dt * dG)
        Y[k+1, IDX_F]     = max(0.0, cF     + dt * dF)
        Y[k+1, IDX_E]     = max(0.0, cE     + dt * dE)
        Y[k+1, IDX_O2]    = max(0.0, cO2    + dt * dO2)
        Y[k+1, IDX_PROT]  = max(0.0, cProt  + dt * dProt)
        Y[k+1, IDX_CARB]  = max(0.0, cCarb  + dt * dCarb)
        Y[k+1, IDX_NREC]  = max(0.0, cNrec  + dt * dNrec)

        for j_aa, ak in enumerate(aa_keys):
            v_aa_net = aa_net_h[ak][k]  # negative = uptake
            if mode_h[-1] == "BIOMASS":
                dAA = cX * v_aa_net  # consumption from medium
            else:
                dAA = 0.0  # extracellular AA unchanged in stationary
            Y[k+1, IDX_AA0 + j_aa] = max(0.0, Y[k, IDX_AA0 + j_aa] + dt * dAA)

        for j_ar, ak_a in enumerate(aroma_keys):
            mw_a = AROMA_MW.get(ak_a, 0.1)
            v_a  = aroma_fh[ak_a][k]
            dA   = mw_a * v_a * cX
            Y[k+1, IDX_AROMA0 + j_ar] = max(0.0, Y[k, IDX_AROMA0 + j_ar] + dt * dA)

    # ══════════════════════════════════════════════════════════════════════════════
    # UNPACK RESULTS
    # ══════════════════════════════════════════════════════════════════════════════

    X_h     = Y[:, IDX_X]
    Nf_h    = Y[:, IDX_NFREE]
    G_h     = Y[:, IDX_G]
    F_h     = Y[:, IDX_F]
    E_h     = Y[:, IDX_E]
    O2_h    = Y[:, IDX_O2]
    Prot_h  = Y[:, IDX_PROT]
    Carb_h  = Y[:, IDX_CARB]
    Nrec_h  = Y[:, IDX_NREC]

    aa_conc_h    = {ak: Y[:, IDX_AA0 + j] for j, ak in enumerate(aa_keys)}
    aroma_conc_h = {ak: Y[:, IDX_AROMA0 + j] for j, ak in enumerate(aroma_keys)}
    aroma_mg_h   = {ak: 1000.0 * aroma_conc_h[ak] for ak in aroma_keys}

    aa_n_h = MW_N * np.sum(
        np.column_stack([aa_conc_h[ak] for ak in aa_keys]), axis=1) if aa_keys else np.zeros(len(tt))
    N_h = Nf_h + aa_n_h

    P_frac_full = Prot_h / np.maximum(X_h, EPS_LOC)
    C_frac_full = Carb_h / np.maximum(X_h, EPS_LOC)

    gam_full = np.array([
        _compute_full_gam(float(P_frac_full[i]), RNA_FRAC, float(C_frac_full[i]),
                          Pbase_global, Rbase_global, Cbase_global)
        for i in range(len(tt))
    ])

    stat_mask = np.array([m == "TURNOVER_ATPM" for m in mode_h], dtype=bool)
    t_stat = float(tt[np.argmax(stat_mask)]) if stat_mask.any() else np.nan

    # ── Mass balance diagnostic ──────────────────────────────────────────────
    # Total N in system = N_free + N_from_AA + N_in_Prot + N_rec
    N_in_prot_h = Y_N_FROM_PROT * Prot_h
    N_total_system_h = Nf_h + aa_n_h + N_in_prot_h + Nrec_h
    N_total_initial  = float(N_total_system_h[0])
    N_total_final    = float(N_total_system_h[-1])
    N_balance_error  = abs(N_total_final - N_total_initial)

    biomass_steps  = sum(m == "BIOMASS" for m in mode_h)
    turnover_steps = sum(m == "TURNOVER_ATPM" for m in mode_h)

    print(f"\n[{scenario_name}]")
    print(f"  Steps: BIOMASS={biomass_steps}, TURNOVER_ATPM={turnover_steps}")
    print(f"  t_stationary ≈ {t_stat:.1f} h")
    print(f"  X_final={X_h[-1]:.4f} gDW/L  |  E_final={E_h[-1]:.4f} g/L")
    print(f"  Prot_final={Prot_h[-1]:.4f} g/L  |  Prot/X={P_frac_full[-1]:.4f}")
    print(f"  N_rec_final={Nrec_h[-1]:.6f} gN/L")
    print(f"  GAM: initial={gam_full[0]:.3f} → final={gam_full[-1]:.3f} mmol ATP/gDW")
    print(f"  N mass balance: initial={N_total_initial:.6f} → final={N_total_final:.6f} gN/L")
    print(f"  N balance error = {N_balance_error:.2e} gN/L  "
          f"({'OK' if N_balance_error < 0.01 else 'CHECK!'})")

    # ══════════════════════════════════════════════════════════════════════════════
    # DataFrames
    # ══════════════════════════════════════════════════════════════════════════════

    states_df = pd.DataFrame({
        "t_h": tt,
        "X_gDW_L": X_h,
        "N_gN_L": N_h,
        "N_free_gN_L": Nf_h,
        "N_from_AA_gN_L": aa_n_h,
        "N_rec_gN_L": Nrec_h,
        "N_in_Prot_gN_L": N_in_prot_h,
        "N_total_system_gN_L": N_total_system_h,
        "G_g_L": G_h,
        "F_g_L": F_h,
        "E_g_L": E_h,
        "O2_g_L": O2_h,
        "Prot_g_L": Prot_h,
        "Carb_g_L": Carb_h,
        "P_frac_gProt_gDW": P_frac_full,
        "C_frac_gCarb_gDW": C_frac_full,
        "GAM_mmol_gDW": gam_full,
    })
    for ak in aa_keys:
        states_df[f"AA_{ak}_mmol_L"] = aa_conc_h[ak]
        states_df[f"AA_{ak}_g_L"]    = aa_conc_h[ak] * AA_MW.get(ak, 0.13)
    for ak in aroma_keys:
        states_df[f"{ak}_g_L"]  = aroma_conc_h[ak]
        states_df[f"{ak}_mg_L"] = aroma_mg_h[ak]

    flux_df = pd.DataFrame({
        "t_h": tt[:-1],
        "mode": mode_h,
        "mu_h_inv": mu_h,
        "v_glu": vglu_h,
        "v_fru": vfru_h,
        "v_eth": veth_h,
        "v_o2": vo2_h,
        "v_N_uptake_gN_gDW_h": vn_h,
        "v_ATPM": vatpm_h,
        "v_Prot": vprot_h,
        "GAM": gam_h,
        "Nrec_inflow_gN_L_h": nrec_inflow_h,
        "Nrec_outflow_gN_L_h": nrec_outflow_h,
        "Nrec_cap_gN_gDW_h": nrec_cap_h,
    })
    for ak in aa_keys:
        flux_df[f"v_uptake_{ak}"] = aa_uptake_h[ak]
    for ak in aroma_keys:
        flux_df[f"v_{ak}"] = aroma_fh[ak]

    display(flux_df.head(10))
    display(flux_df.tail(10))

    # ══════════════════════════════════════════════════════════════════════════════
    # PLOTS
    # ══════════════════════════════════════════════════════════════════════════════

    # ── 1. Main states ──────────────────────────────────────────────────────
    state_series = [
        ("X (gDW/L)",      X_h,     mu_h,    "mu (1/h)"),
        ("N_ext (gN/L)",   N_h,     vn_h,    "v_N"),
        ("N_rec (gN/L)",   Nrec_h,  nrec_inflow_h, "Nrec inflow"),
        ("Glucosa (g/L)",  G_h,     vglu_h,  "v_glu"),
        ("Fructosa (g/L)", F_h,     vfru_h,  "v_fru"),
        ("Etanol (g/L)",   E_h,     veth_h,  "v_eth"),
        ("O2 (g/L)",       O2_h,    vo2_h,   "v_o2"),
        ("Prot (g/L)",     Prot_h,  vprot_h, "v_Prot"),
        ("Carb (g/L)",     Carb_h,  np.zeros(n_st), "—"),
        ("N total sistema (gN/L)", N_total_system_h, np.zeros(n_st), "—"),
    ] + [
        (f"AA_{ak} (mmol/L)", aa_conc_h[ak], aa_uptake_h[ak], f"v_{ak}")
        for ak in aa_keys
    ] + [
        (f"{ak} (mg/L)", aroma_mg_h[ak], aroma_fh[ak], f"v_{ak}")
        for ak in aroma_keys
    ]

    ncols_s = 4
    nrows_s = math.ceil(len(state_series) / ncols_s)
    fig, ax = plt.subplots(nrows_s, ncols_s, figsize=(5.2*ncols_s, 3.5*nrows_s),
                           sharex=True, dpi=120)
    ax = np.atleast_1d(ax).ravel()
    for i, (ttl, yy_, vv_, vlbl) in enumerate(state_series):
        ax[i].plot(tt, yy_, lw=2, color="tab:blue")
        if np.isfinite(t_stat):
            ax[i].axvline(t_stat, ls="--", lw=1, color="tab:blue", alpha=0.5)
        ax[i].set_title(ttl, fontsize=9)
        ax[i].grid(alpha=0.3)
        ax[i].set_xlabel("t (h)", fontsize=8)
        ax2 = ax[i].twinx()
        ax2.step(tt[:-1], vv_, where="post", lw=1.4, ls="--", color="tab:orange", alpha=0.85)
        ax2.set_ylabel(vlbl, color="tab:orange", fontsize=7)
        ax2.tick_params(axis="y", labelcolor="tab:orange", labelsize=7)
    for j in range(i+1, len(ax)):
        ax[j].axis("off")
    plt.suptitle(f"Estados dFBA – {scenario_name}", fontsize=11, y=1.01)
    plt.tight_layout()
    plt.show()

    # ── 2. GAM + composition ──────────────────────────────────────────────────
    fig_gam, ax_gam = plt.subplots(1, 3, figsize=(15, 4), dpi=120)

    ax_gam[0].plot(tt, gam_full, lw=2.5, color="tab:red")
    ax_gam[0].axhline(30.49, ls="--", lw=1.2, color="gray", label="GAM_base")
    if np.isfinite(t_stat):
        ax_gam[0].axvline(t_stat, ls=":", lw=1.5, color="navy", alpha=0.6)
    ax_gam[0].set_title("GAM (mmol ATP/gDW)")
    ax_gam[0].set_xlabel("t (h)")
    ax_gam[0].legend(fontsize=8)
    ax_gam[0].grid(alpha=0.3)

    ax_gam[1].plot(tt, P_frac_full, lw=2.5, color="tab:green", label="Protein (g/gDW)")
    ax_gam[1].plot(tt, C_frac_full, lw=2.5, color="tab:purple", label="Carb (g/gDW)")
    ax_gam[1].axhline(Pbase_global, ls="--", lw=1, color="tab:green", alpha=0.5)
    ax_gam[1].axhline(Cbase_global, ls="--", lw=1, color="tab:purple", alpha=0.5)
    if np.isfinite(t_stat):
        ax_gam[1].axvline(t_stat, ls=":", lw=1.5, color="navy", alpha=0.6)
    ax_gam[1].set_title("Composición celular (g/gDW)")
    ax_gam[1].set_xlabel("t (h)")
    ax_gam[1].legend(fontsize=8)
    ax_gam[1].grid(alpha=0.3)

    ax_gam[2].plot(tt, Prot_h, lw=2.5, color="tab:green", label="Prot (g/L)")
    ax_gam[2].plot(tt, X_h,    lw=2.5, color="tab:blue",  label="X (gDW/L)")
    ax_gam[2].plot(tt, Nrec_h, lw=2.5, color="tab:red",   label="N_rec (gN/L)")
    if np.isfinite(t_stat):
        ax_gam[2].axvline(t_stat, ls=":", lw=1.5, color="navy", alpha=0.6)
    ax_gam[2].set_title("Biomasa, Proteína y N_rec")
    ax_gam[2].set_xlabel("t (h)")
    ax_gam[2].legend(fontsize=8)
    ax_gam[2].grid(alpha=0.3)

    plt.suptitle(f"GAM variable y N_rec – {scenario_name}", fontsize=11)
    plt.tight_layout()
    plt.show()

    # ── 3. N_rec inflow/outflow ──────────────────────────────────────────────
    fig_nr, ax_nr = plt.subplots(1, 3, figsize=(15, 4), dpi=120)

    ax_nr[0].plot(tt, Nrec_h, lw=2.5, color="tab:red")
    ax_nr[0].set_title("N_rec pool (gN/L)")
    ax_nr[0].set_xlabel("t (h)")
    ax_nr[0].grid(alpha=0.3)
    if np.isfinite(t_stat):
        ax_nr[0].axvline(t_stat, ls=":", lw=1.5, color="navy", alpha=0.6)

    ax_nr[1].step(tt[:-1], nrec_inflow_h, where="post", lw=2, color="tab:green", label="inflow (from Prot)")
    ax_nr[1].step(tt[:-1], nrec_outflow_h, where="post", lw=2, ls="--", color="tab:orange", label="outflow (to FBA)")
    ax_nr[1].set_title("N_rec flows (gN/L/h)")
    ax_nr[1].set_xlabel("t (h)")
    ax_nr[1].legend(fontsize=8)
    ax_nr[1].grid(alpha=0.3)
    if np.isfinite(t_stat):
        ax_nr[1].axvline(t_stat, ls=":", lw=1.5, color="navy", alpha=0.6)

    ax_nr[2].plot(tt, N_total_system_h, lw=2.5, color="tab:brown")
    ax_nr[2].set_title("N total sistema (gN/L) – balance check")
    ax_nr[2].set_xlabel("t (h)")
    ax_nr[2].grid(alpha=0.3)
    y_range = max(0.001, float(np.ptp(N_total_system_h)))
    ax_nr[2].set_ylim(N_total_initial - 2*y_range, N_total_initial + 2*y_range)

    plt.suptitle(f"Balance de N_rec – {scenario_name}", fontsize=11)
    plt.tight_layout()
    plt.show()

    # ── 4. Aromas ────────────────────────────────────────────────────────────
    if aroma_keys:
        ncols_a = min(len(aroma_keys), 5)
        nrows_a = math.ceil(len(aroma_keys) / ncols_a)
        fig_a, ax_a = plt.subplots(nrows_a, ncols_a,
                                    figsize=(4.5*ncols_a, 3.5*nrows_a),
                                    dpi=120, squeeze=False)
        ax_a = ax_a.ravel()
        for i_a, ak_a in enumerate(aroma_keys):
            ax_a[i_a].plot(tt, aroma_mg_h[ak_a], lw=2.5, color="tab:brown")
            ax_a[i_a].set_title(f"{ak_a} (mg/L)", fontsize=9)
            ax_a[i_a].set_xlabel("t (h)")
            ax_a[i_a].grid(alpha=0.3)
            if np.isfinite(t_stat):
                ax_a[i_a].axvline(t_stat, ls=":", lw=1.2, color="navy", alpha=0.5)
        for j_a in range(i_a+1, len(ax_a)):
            ax_a[j_a].axis("off")
        plt.suptitle(f"Aromas – {scenario_name}", fontsize=11)
        plt.tight_layout()
        plt.show()

    # ── Save CSV ─────────────────────────────────────────────────────────────
    sfx_file = "aerobic_nrec" if aerobic_mode else "anaerobic_nrec"
    out_states = OUT_DIR / f"dfba_nrec_states_{sfx_file}.csv"
    out_flux   = OUT_DIR / f"dfba_nrec_flux_{sfx_file}.csv"
    states_df.to_csv(out_states, index=False)
    flux_df.to_csv(out_flux, index=False)
    print(f"\nArchivos guardados:\n  {out_states}\n  {out_flux}")

    return {
        "states_df":       states_df,
        "flux_df":         flux_df,
        "Y":               Y,
        "mode_h":          mode_h,
        "gam_full":        gam_full,
        "P_frac_full":     P_frac_full,
        "C_frac_full":     C_frac_full,
        "time_grid_nrec":  tt,
        "N_total_system":  N_total_system_h,
        "N_balance_error": N_balance_error,
    }


# ══════════════════════════════════════════════════════════════════════════════════════
# RUN
# ══════════════════════════════════════════════════════════════════════════════════════

res_nrec = _run_dfba_nrec(
    scenario_name="Anaeróbico + N_rec conservativo",
    aerobic_mode=True,
    o2_init_g_l=0.008,
    dfba_dt=1.0,
    t_end=72.0,
)

# ── Update globals so downstream cells can use the new trajectory ──────────────────
# Overwrite states_df_phase and flux_df_phase so that:
#   - _build_phase_snapshot_model works on the new trajectory
#   - acetate ester analysis cells pick up the new states
states_df_phase = res_nrec["states_df"].copy()
flux_df_phase   = res_nrec["flux_df"].copy()

print("\n[INFO] states_df_phase y flux_df_phase actualizados con trayectoria N_rec.")
print("       Las celdas downstream de snapshots/FVA usarán esta trayectoria.")