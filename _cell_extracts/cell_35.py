# ======================================================================================
# dFBA SOA + pFBA con protein-turnover y tracking de aromas fermentativos
# ======================================================================================

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cobra

print("\n" + "=" * 100)
print("dFBA SOA + pFBA con protein-turnover + aromas")
print("=" * 100)

required_symbols = [
    "model", "rxn_by_id", "required_ids", "base_lb", "base_ub",
    "KINETIC_N_SOURCE_IDS", "N_VITAMIN_EX_IDS", "N_atoms_map",
    "OBJ_ID", "GLU_ID", "FRU_ID", "ETH_ID", "O2_ID", "ATPM_ID",
    "ATPM_LB_NO_GROWTH", "ATPM_UB_NO_GROWTH",
    "_kinetic_limits", "DFBA_DT", "time_grid",
    "N_DEPLETION_CLOSE_THRESHOLD", "O2_DEPLETION_THRESHOLD",
    "O2_VMAX_UPTAKE", "O2_KO_G_L",
    "MW_GLU", "MW_FRU", "MW_ETH", "MW_N", "MW_O2",
    "OUT_DIR", "APPLY_PRODUCT_CAPS", "EPS",
]
missing_symbols = [s for s in required_symbols if s not in globals()]
if missing_symbols:
    raise RuntimeError(f"Faltan símbolos requeridos: {missing_symbols}")


def _rxn_search_text(rxn):
    met_names = " ".join(getattr(m, "name", "") for m in rxn.metabolites)
    return f"{rxn.id} {rxn.name} {met_names}".lower()


def _is_exchange_like(rxn):
    nm = (rxn.name or "").lower()
    return ("exchange" in nm) or (len(rxn.metabolites) == 1)


def _find_candidates(keyword_sets, exchange_preferred=True):
    found = []
    for rxn in model.reactions:
        txt = _rxn_search_text(rxn)
        for kws in keyword_sets:
            if all(str(k).lower() in txt for k in kws):
                score = 0
                if exchange_preferred and _is_exchange_like(rxn):
                    score -= 10
                score += len(rxn.metabolites)
                found.append((score, rxn.id, rxn.name))
                break
    return sorted(set(found))


def _pick_first(keyword_sets, manual_id=None, exchange_preferred=True):
    if manual_id is not None:
        if manual_id not in {r.id for r in model.reactions}:
            raise RuntimeError(f"Manual ID no encontrado: {manual_id}")
        return manual_id, [(0, manual_id, model.reactions.get_by_id(manual_id).name)]
    cands = _find_candidates(keyword_sets, exchange_preferred=exchange_preferred)
    rid = cands[0][1] if cands else None
    return rid, cands


MANUAL_PROT_RXN_ID = "r_4047"
MANUAL_AA_EX_IDS = {"phe": None, "leu": None, "val": None, "met": None, "tyr": None}
MANUAL_AROMA_EX_IDS = {"pea": None, "isoamyl": None, "isobutanol": None, "methionol": None, "tyrosol": None, "ethyl_acetate": None}

PROT_SEARCH = [["protein", "exchange"], ["protein", "production"], ["protein"]]
AA_SEARCH = {
    "phe": [["phenylalanine", "exchange"], ["l-phenylalanine"], ["phenylalanine"]],
    "leu": [["leucine", "exchange"], ["l-leucine"], ["leucine"]],
    "val": [["valine", "exchange"], ["l-valine"], ["valine"]],
    "met": [["methionine", "exchange"], ["l-methionine"], ["methionine"]],
    "tyr": [["tyrosine", "exchange"], ["l-tyrosine"], ["tyrosine"]],
}
AROMA_SEARCH = {
    "pea": [["2-phenylethanol", "exchange"], ["phenylethanol", "exchange"], ["pea", "exchange"]],
    "isoamyl": [["isoamyl alcohol", "exchange"], ["3-methyl-1-butanol", "exchange"], ["isoamyl"]],
    "isobutanol": [["isobutanol", "exchange"], ["2-methyl-1-propanol", "exchange"], ["isobutanol"]],
    "methionol": [["methionol", "exchange"], ["3-methylthio-1-propanol", "exchange"], ["methionol"]],
    "tyrosol": [["tyrosol", "exchange"], ["tyrosol"]],
    "ethyl_acetate": [["ethyl acetate", "exchange"], ["ethyl-acetate", "exchange"], ["ethyl acetate"]],
}

PROT_RXN_ID, prot_cands = _pick_first(PROT_SEARCH, manual_id=MANUAL_PROT_RXN_ID, exchange_preferred=True)
AA_EX_IDS, AROMA_EX_IDS = {}, {}
for k, qs in AA_SEARCH.items():
    AA_EX_IDS[k], _ = _pick_first(qs, manual_id=MANUAL_AA_EX_IDS.get(k), exchange_preferred=True)
for k, qs in AROMA_SEARCH.items():
    AROMA_EX_IDS[k], _ = _pick_first(qs, manual_id=MANUAL_AROMA_EX_IDS.get(k), exchange_preferred=True)

found_precursors = {k: v for k, v in AA_EX_IDS.items() if v is not None}
found_aromas = {k: v for k, v in AROMA_EX_IDS.items() if v is not None}
if PROT_RXN_ID is None or len(found_precursors) < 3:
    raise RuntimeError("Revisa mapeo PROT/AA en MANUAL_*_IDS")

# ---------------------------
# Parámetros
# ---------------------------
PROT_CONTENT = 0.45
K_DEATH = 0.005
TURNOVER_LAMBDA = 0.03
VPROT_FLOOR_FRAC = 0.90
XA_FRACTION = 1.0
N_AMMONIA_FRACTION = 0.50  # fracción de N macroscópico inicial como N amoniacal
N_TOTAL_DEPLETION_THRESHOLD = 1e-3

# NUEVO: cinética de uptake en crecimiento (misma k para todos)
K_AA_UPTAKE_GROWTH = 0.08  # 1/h
BIOMASS_LOCK_FRAC = 0.999  # lock de biomasa para fase B en crecimiento
AA_UPTAKE_WEIGHT = 1.0

# Refinamiento temporal del dFBA
DFBA_DT_FINE = min(float(DFBA_DT), 1.0)
TIME_GRID_FINE = np.arange(float(time_grid[0]), float(time_grid[-1]) + 0.5 * DFBA_DT_FINE, DFBA_DT_FINE)

AA_ALPHA_TOTAL = 1.0
AA_ALPHA_STATIC = AA_ALPHA_TOTAL / max(1, len(found_precursors))
AA_ALPHA = {k: AA_ALPHA_STATIC for k in found_precursors.keys()}

AA_MW = {
    "phe": 0.16519,
    "leu": 0.13117,
    "val": 0.11715,
    "met": 0.14921,
    "tyr": 0.18119,
}
AROMA_MW = {
    "pea": 0.12217,
    "isoamyl": 0.08815,
    "isobutanol": 0.07412,
    "methionol": 0.10619,
    "tyrosol": 0.13816,
    "ethyl_acetate": 0.08811,
}

TRACKED_RXN_IDS = list(dict.fromkeys(
    list(required_ids)
    + [PROT_RXN_ID]
    + [rid for rid in AA_EX_IDS.values() if rid is not None]
    + [rid for rid in AROMA_EX_IDS.values() if rid is not None]
) )
base_lb_turn, base_ub_turn = {}, {}
for rid in TRACKED_RXN_IDS:
    rxn = model.reactions.get_by_id(rid)
    base_lb_turn[rid] = float(rxn.lower_bound)
    base_ub_turn[rid] = float(rxn.upper_bound)


def _run_dfba_turnover_aroma(*, scenario_name="Anaeróbico-turnover", aerobic_mode=False, o2_init_g_l=0.0, dfba_dt=DFBA_DT_FINE, time_grid_local=TIME_GRID_FINE):
    aroma_keys = list(found_aromas.keys())
    aa_keys = list(found_precursors.keys())
    n_aroma = len(aroma_keys)
    n_aa = len(aa_keys)
    dt = float(dfba_dt)
    tt = np.asarray(time_grid_local, dtype=float)

    # Partición de N inicial: 50% amoniacal + 50% en AA (uniforme entre AA_i)
    N0_total = float(globals().get("N0", 0.14))
    N0_ammonia = N0_total * N_AMMONIA_FRACTION
    N0_from_aa = max(0.0, N0_total - N0_ammonia)
    aa0_each_mmol_L = (N0_from_aa / MW_N) / max(1, n_aa)

    # Estados: [X,N,G,F,E,O2,Prot,AA_i...,aromas...]
    Y = np.zeros((len(tt), 7 + n_aa + n_aroma), dtype=float)
    Y[0, 0] = 0.5
    Y[0, 1] = N0_ammonia
    Y[0, 2] = 110.0
    Y[0, 3] = 110.0
    Y[0, 4] = 0.0
    Y[0, 5] = o2_init_g_l
    Y[0, 6] = Y[0, 0] * PROT_CONTENT # Se asueme esta composición proteica fija de la biomasa, para simplificar el tracking de turnover proteico
    for idx in range(n_aa):
        Y[0, 7 + idx] = aa0_each_mmol_L #Creamos estados de modelo para cada AA precursor, con la concentración inicial correspondiente a la fracción de N asignada a AA

    n_steps = len(tt) - 1
    mu_hist = np.zeros(n_steps)
    vglu_hist, vfru_hist, veth_hist = np.zeros(n_steps), np.zeros(n_steps), np.zeros(n_steps)
    vo2_hist, vn_eff_hist = np.zeros(n_steps), np.zeros(n_steps)
    vatpm_hist, vprot_hist, vprot_star_hist = np.zeros(n_steps), np.zeros(n_steps), np.zeros(n_steps)
    turn_cap_hist, turn_uptake_hist = np.zeros(n_steps), np.zeros(n_steps)
    mode_hist = []

    aroma_flux_hist = {k: np.zeros(n_steps) for k in aroma_keys}
    precursor_uptake_hist = {k: np.zeros(n_steps) for k in aa_keys}
    precursor_net_flux_hist = {k: np.zeros(n_steps) for k in aa_keys}
    aa_release_hist = {k: np.zeros(n_steps) for k in aa_keys}
    aa_turnover_internal_hist = {k: np.zeros(n_steps) for k in aa_keys}

    interest_rxn_ids = [OBJ_ID, GLU_ID, FRU_ID, ETH_ID, O2_ID, ATPM_ID, PROT_RXN_ID] + list(found_precursors.values()) + list(found_aromas.values())
    interest_rxn_ids = list(dict.fromkeys(interest_rxn_ids))
    v_interest_hist = {rid: np.zeros(n_steps) for rid in interest_rxn_ids}

    def _apply_bounds_turnover(mtmp, y):
        cX, cN_free, cG, cF, cE, cO2, cProt = [max(0.0, float(v)) for v in y[:7]]
        aa_state = {aa_key: max(0.0, float(y[7 + idx])) for idx, aa_key in enumerate(aa_keys)}
        aa_n_total = MW_N * float(sum(aa_state.values()))
        cN_total = cN_free + aa_n_total

        lim_glu, lim_fru, lim_eth, lim_obj, lim_n_map = _kinetic_limits(cX, cN_free, cG, cF, cE, cO2) #obtenemos límites cinéticos dinámicos de modelo zenteno
        rb = {r.id: r for r in mtmp.reactions}

        for rid in TRACKED_RXN_IDS:
            rb[rid].lower_bound = base_lb_turn[rid]
            rb[rid].upper_bound = base_ub_turn[rid]

        rb[GLU_ID].lower_bound = max(rb[GLU_ID].lower_bound, -lim_glu)
        rb[GLU_ID].upper_bound = min(rb[GLU_ID].upper_bound, 0.0)
        rb[FRU_ID].lower_bound = max(rb[FRU_ID].lower_bound, -lim_fru)
        rb[FRU_ID].upper_bound = min(rb[FRU_ID].upper_bound, 0.0)

        for rid in KINETIC_N_SOURCE_IDS:
            Li = float(lim_n_map[rid])
            rb[rid].lower_bound = max(rb[rid].lower_bound, -Li)
            rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)

        if APPLY_PRODUCT_CAPS:
            rb[ETH_ID].upper_bound = min(rb[ETH_ID].upper_bound, lim_eth)
            rb[OBJ_ID].upper_bound = min(rb[OBJ_ID].upper_bound, lim_obj)
        else:
            rb[ETH_ID].upper_bound = max(rb[ETH_ID].upper_bound, 1000.0)
            rb[OBJ_ID].upper_bound = max(rb[OBJ_ID].upper_bound, 1000.0)

        if aerobic_mode:
            lim_o2 = O2_VMAX_UPTAKE * (cO2 / (cO2 + O2_KO_G_L + EPS))
            lim_o2 = max(0.0, float(lim_o2))
            if cO2 <= O2_DEPLETION_THRESHOLD:
                lim_o2 = 0.0
            rb[O2_ID].lower_bound = -lim_o2
            rb[O2_ID].upper_bound = 0.0
        else:
            rb[O2_ID].lower_bound = 0.0
            rb[O2_ID].upper_bound = 0.0

        turnover_caps = {k: 0.0 for k in aa_keys}

        if cN_total > N_TOTAL_DEPLETION_THRESHOLD:
            for aa_key, rid in found_precursors.items():
                q_i = K_AA_UPTAKE_GROWTH * aa_state[aa_key] / max(cX, EPS)
                q_i = max(0.0, float(q_i))
                rb[rid].lower_bound = -q_i
                rb[rid].upper_bound = min(float(rb[rid].upper_bound), 0.0)
                turnover_caps[aa_key] = q_i
            qturn_total = float(sum(turnover_caps.values()))
        else:
            for rid in KINETIC_N_SOURCE_IDS:
                rb[rid].lower_bound = max(float(rb[rid].lower_bound), 0.0)
                rb[rid].upper_bound = min(float(rb[rid].upper_bound), 0.0)
            for rid in N_VITAMIN_EX_IDS:
                if rid in rb:
                    rb[rid].lower_bound = max(float(rb[rid].lower_bound), 0.0)
                    rb[rid].upper_bound = min(float(rb[rid].upper_bound), 0.0)

            XA = max(EPS, XA_FRACTION * cX)
            for aa_key, rid in found_precursors.items():
                q_i = TURNOVER_LAMBDA * cProt * AA_ALPHA[aa_key] / XA
                q_i = max(0.0, float(q_i))
                rb[rid].lower_bound = -q_i
                rb[rid].upper_bound = min(float(rb[rid].upper_bound), 0.0)
                turnover_caps[aa_key] = q_i
            qturn_total = float(sum(turnover_caps.values()))

        return rb, turnover_caps, qturn_total, cN_total

    for k in range(n_steps):
        yk = np.maximum(Y[k, :], 0.0)
        cX, cN, cG, cF, cE, cO2, cProt = yk[:7]

        with model as mtmp:
            rb, turnover_caps, qturn_total, cN_total = _apply_bounds_turnover(mtmp, yk)

            if cN_total > N_TOTAL_DEPLETION_THRESHOLD:
                mode_hist.append("BIOMASS")
                try:
                    mtmp.objective = OBJ_ID
                    mtmp.objective_direction = "max"
                    sol_mu = mtmp.optimize()
                except Exception:
                    mode_hist[-1] = "FAIL"
                    continue
                if getattr(sol_mu, "status", "") != "optimal":
                    mode_hist[-1] = "NO_OPT"
                    continue

                mu_star = max(0.0, float(sol_mu.fluxes.get(OBJ_ID, 0.0)))
                rb[OBJ_ID].lower_bound = max(float(rb[OBJ_ID].lower_bound), BIOMASS_LOCK_FRAC * mu_star)

                uptake_obj = {mtmp.reactions.get_by_id(rid): -AA_UPTAKE_WEIGHT for rid in found_precursors.values()}
                mtmp.objective = uptake_obj
                mtmp.objective_direction = "max"
                try:
                    sol = cobra.flux_analysis.pfba(mtmp, fraction_of_optimum=1.0)
                except Exception:
                    mode_hist[-1] = "FAIL"
                    continue
                vprot_star = 0.0
            else:
                rb[OBJ_ID].lower_bound = 0.0
                rb[OBJ_ID].upper_bound = 0.0

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
                    rb[PROT_RXN_ID].lower_bound = max(float(rb[PROT_RXN_ID].lower_bound), VPROT_FLOOR_FRAC * vprot_star)

                rb[ATPM_ID].lower_bound = ATPM_LB_NO_GROWTH
                rb[ATPM_ID].upper_bound = ATPM_UB_NO_GROWTH

                mtmp.objective = ATPM_ID
                mtmp.objective_direction = "max"
                mode_hist.append("TURNOVER_ATPM")

                try:
                    sol = cobra.flux_analysis.pfba(mtmp, fraction_of_optimum=1.0)
                except Exception:
                    mode_hist[-1] = "FAIL"
                    continue

            if getattr(sol, "status", "optimal") != "optimal":
                mode_hist[-1] = "NO_OPT"
                continue

            v_obj = float(sol.fluxes.get(OBJ_ID, 0.0))
            v_glu = float(sol.fluxes.get(GLU_ID, 0.0))
            v_fru = float(sol.fluxes.get(FRU_ID, 0.0))
            v_eth = float(sol.fluxes.get(ETH_ID, 0.0))
            v_o2 = float(sol.fluxes.get(O2_ID, 0.0))
            v_atpm = float(sol.fluxes.get(ATPM_ID, 0.0))
            v_prot = float(sol.fluxes.get(PROT_RXN_ID, 0.0))

            if cN_total > N_TOTAL_DEPLETION_THRESHOLD:
                vn_eff = 0.0
                for rid in KINETIC_N_SOURCE_IDS:
                    v_i = float(sol.fluxes.get(rid, 0.0))
                    vn_eff += max(0.0, -v_i) * N_atoms_map[rid] * MW_N
            else:
                vn_eff = 0.0

            turn_uptake = 0.0
            for aa_key, rid in found_precursors.items():
                v_i = float(sol.fluxes.get(rid, 0.0))
                q_i = max(0.0, -v_i)
                precursor_uptake_hist[aa_key][k] = q_i
                precursor_net_flux_hist[aa_key][k] = v_i
                aa_turnover_internal_hist[aa_key][k] = TURNOVER_LAMBDA * cProt * AA_ALPHA[aa_key] if cN_total <= N_TOTAL_DEPLETION_THRESHOLD else 0.0
                aa_release_hist[aa_key][k] = 0.0
                turn_uptake += q_i

            for akey in aroma_keys:
                rid = found_aromas[akey]
                v_a = float(sol.fluxes.get(rid, 0.0))
                aroma_flux_hist[akey][k] = max(0.0, v_a)

            for rid in interest_rxn_ids:
                v_interest_hist[rid][k] = float(sol.fluxes.get(rid, 0.0))

        mu_hist[k] = max(0.0, v_obj)
        vglu_hist[k] = max(0.0, -v_glu)
        vfru_hist[k] = max(0.0, -v_fru)
        veth_hist[k] = max(0.0, v_eth)
        vo2_hist[k] = max(0.0, -v_o2)
        vn_eff_hist[k] = max(0.0, vn_eff)
        vatpm_hist[k] = max(0.0, v_atpm)
        vprot_hist[k] = max(0.0, v_prot)
        vprot_star_hist[k] = max(0.0, vprot_star)
        turn_cap_hist[k] = max(0.0, qturn_total)
        turn_uptake_hist[k] = max(0.0, turn_uptake)

        dX = v_obj * cX if mode_hist[-1] == "BIOMASS" else 0.0
        dN = -vn_eff * cX if mode_hist[-1] == "BIOMASS" else 0.0
        dG = -MW_GLU * max(0.0, -v_glu) * cX
        dF = -MW_FRU * max(0.0, -v_fru) * cX
        dE = MW_ETH * max(0.0, v_eth) * cX
        dO2 = -MW_O2 * max(0.0, -v_o2) * cX

        Xv = cX
        Xa = XA_FRACTION * cX
        dProt = (
            max(0.0, v_obj) * Xv * PROT_CONTENT
            - cProt * K_DEATH
            + Xa * max(0.0, v_prot)
            - TURNOVER_LAMBDA * cProt
        )

        Y[k + 1, 0] = max(0.0, cX + dt * dX)
        Y[k + 1, 1] = max(0.0, cN + dt * dN)
        Y[k + 1, 2] = max(0.0, cG + dt * dG)
        Y[k + 1, 3] = max(0.0, cF + dt * dF)
        Y[k + 1, 4] = max(0.0, cE + dt * dE)
        Y[k + 1, 5] = max(0.0, cO2 + dt * dO2)
        Y[k + 1, 6] = max(0.0, cProt + dt * dProt)

        aa_start = 7
        for j, aa_key in enumerate(aa_keys, start=aa_start):
            v_aa_i = precursor_net_flux_hist[aa_key][k]
            rel_aa_i = aa_release_hist[aa_key][k]
            dAA_i = cX * v_aa_i + rel_aa_i if mode_hist[-1] == "BIOMASS" else 0.0
            Y[k + 1, j] = max(0.0, Y[k, j] + dt * dAA_i)

        aroma_start = 7 + n_aa
        for j, akey in enumerate(aroma_keys, start=aroma_start):
            mw = AROMA_MW[akey]
            v_a = aroma_flux_hist[akey][k]
            dA = mw * v_a * cX
            Y[k + 1, j] = max(0.0, Y[k, j] + dt * dA)

    X_hist, N_free_hist, G_hist, F_hist, E_hist, O2_hist, Prot_hist = Y[:, 0], Y[:, 1], Y[:, 2], Y[:, 3], Y[:, 4], Y[:, 5], Y[:, 6]
    aa_conc_hist = {aa_key: Y[:, 7 + idx] for idx, aa_key in enumerate(aa_keys)}
    aa_conc_hist_g = {aa_key: aa_conc_hist[aa_key] * AA_MW[aa_key] for aa_key in aa_keys}
    aa_n_hist = MW_N * np.sum(np.column_stack([aa_conc_hist[aa_key] for aa_key in aa_keys]), axis=1)
    N_hist = N_free_hist + aa_n_hist
    aroma_conc_hist = {akey: Y[:, 7 + n_aa + j] for j, akey in enumerate(aroma_keys)}
    aroma_conc_hist_mg = {akey: 1000.0 * aroma_conc_hist[akey] for akey in aroma_keys}

    stationary_mask = np.array([m == "TURNOVER_ATPM" for m in mode_hist], dtype=bool)
    t_stationary = float(tt[np.argmax(stationary_mask)]) if stationary_mask.any() else np.nan

    states_df = pd.DataFrame({
        "t_h": tt,
        "X_gDW_L": X_hist,
        "N_gN_L": N_hist,
        "N_free_gN_L": N_free_hist,
        "N_from_AA_gN_L": aa_n_hist,
        "G_g_L": G_hist,
        "F_g_L": F_hist,
        "E_g_L": E_hist,
        "O2_g_L": O2_hist,
        "Prot_g_L": Prot_hist,
    })
    for aa_key in aa_keys:
        states_df[f"AA_{aa_key}_mmol_L"] = aa_conc_hist[aa_key]
        states_df[f"AA_{aa_key}_g_L"] = aa_conc_hist_g[aa_key]
    for akey in aroma_keys:
        states_df[f"{akey}_g_L"] = aroma_conc_hist[akey]
        states_df[f"{akey}_mg_L"] = aroma_conc_hist_mg[akey]

    flux_df = pd.DataFrame({
        "t_h": tt[:-1],
        "mode": mode_hist,
        "mu_h_inv": mu_hist,
        "v_glu_uptake_mmol_gDW_h": vglu_hist,
        "v_fru_uptake_mmol_gDW_h": vfru_hist,
        "v_eth_prod_mmol_gDW_h": veth_hist,
        "v_o2_uptake_mmol_gDW_h": vo2_hist,
        "v_N_uptake_gN_gDW_h": vn_eff_hist,
        "v_ATPM_mmol_gDW_h": vatpm_hist,
        "v_Prot_mmol_gDW_h": vprot_hist,
        "v_Prot_star_mmol_gDW_h": vprot_star_hist,
        "turnover_cap_total_mmol_gDW_h": turn_cap_hist,
        "turnover_uptake_total_mmol_gDW_h": turn_uptake_hist,
    })
    for aa_key in aa_keys:
        flux_df[f"v_uptake_{aa_key}_mmol_gDW_h"] = precursor_uptake_hist[aa_key]
        flux_df[f"v_release_{aa_key}_mmol_L_h"] = aa_release_hist[aa_key]
        flux_df[f"v_turnover_internal_{aa_key}_mmol_L_h"] = aa_turnover_internal_hist[aa_key]
        flux_df[f"v_turnover_{aa_key}_mmol_gDW_h"] = aa_turnover_internal_hist[aa_key] / np.maximum(X_hist[:-1], EPS)
    for akey in aroma_keys:
        flux_df[f"v_{akey}_prod_mmol_gDW_h"] = aroma_flux_hist[akey]

    rxn_labels = {
        OBJ_ID: "mu",
        GLU_ID: "v_glu",
        FRU_ID: "v_fru",
        ETH_ID: "v_eth",
        O2_ID: "v_o2",
        ATPM_ID: "v_ATPM",
        PROT_RXN_ID: "v_Prot",
    }
    for aa_key, rid in found_precursors.items():
        rxn_labels[rid] = f"v_{aa_key}"
    for akey, rid in found_aromas.items():
        rxn_labels[rid] = f"v_{akey}"

    v_interest_df = pd.DataFrame({"t_h": tt[:-1], "mode": mode_hist})
    for rid in interest_rxn_ids:
        v_interest_df[f"{rxn_labels.get(rid, rid)}_mmol_gDW_h"] = v_interest_hist[rid]

    print("\n[Indicadores globales]")
    display(pd.DataFrame([{
        "scenario": scenario_name,
        "dfba_dt_h": dt,
        "n_dfba_steps": n_steps,
        "t_stationary_h": t_stationary,
        "X_final_gDW_L": float(X_hist[-1]),
        "Prot_final_g_L": float(Prot_hist[-1]),
        "N0_total_gN_L": N0_total,
        "N0_ammonia_gN_L": N0_ammonia,
        "AA0_each_mmol_L": aa0_each_mmol_L,
        "N_total_depletion_threshold_gN_L": N_TOTAL_DEPLETION_THRESHOLD,
        "K_AA_uptake_growth_h_inv": K_AA_UPTAKE_GROWTH,
        "AA_alpha_static": AA_ALPHA_STATIC,
        "AA_uptake_growth_mean_mmol_gDW_h": float(np.mean(turn_uptake_hist[np.array([m == "BIOMASS" for m in mode_hist])])) if any(m == "BIOMASS" for m in mode_hist) else np.nan,
        "turnover_uptake_mean_stat": float(np.mean(turn_uptake_hist[stationary_mask])) if stationary_mask.any() else np.nan,
    }]))
    print("\n[Flujos relevantes en el tiempo]")
    display(flux_df)
    print("\n[Solución de v en el tiempo]")
    display(v_interest_df)

    state_series = [
        ("X (gDW/L)", X_hist, mu_hist, "mu (1/h)"),
        ("N (gN/L)", N_hist, vn_eff_hist, "v_N (gN/gDW/h)"),
        ("Glucosa (g/L)", G_hist, vglu_hist, "v_glu (mmol/gDW/h)"),
        ("Fructosa (g/L)", F_hist, vfru_hist, "v_fru (mmol/gDW/h)"),
        ("Etanol (g/L)", E_hist, veth_hist, "v_eth (mmol/gDW/h)"),
        ("O2 (g/L)", O2_hist, vo2_hist, "v_o2 (mmol/gDW/h)"),
        ("Prot (g/L)", Prot_hist, vprot_hist, "v_Prot (mmol/gDW/h)"),
    ] + [
        (f"AA_{aa_key} (mmol/L)", aa_conc_hist[aa_key], precursor_uptake_hist[aa_key], f"v_{aa_key} (mmol/gDW/h)")
        for aa_key in aa_keys
    ] + [
        (f"{akey} (mg/L)", aroma_conc_hist_mg[akey], aroma_flux_hist[akey], f"v_{akey} (mmol/gDW/h)")
        for akey in aroma_keys
    ]

    ncols, nrows = 4, math.ceil(len(state_series) / 4)
    fig, ax = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 3.6 * nrows), sharex=True, dpi=160)
    ax = np.atleast_1d(ax).ravel()
    for i, (ttl, yy, vv, vlabel) in enumerate(state_series):
        ax[i].plot(tt, yy, lw=2, color="tab:blue")
        if np.isfinite(t_stationary):
            ax[i].axvline(t_stationary, ls="--", lw=1, color="tab:blue")
        if ttl.startswith("AA_"):
            ymax = max(0.05, 1.10 * float(np.max(yy)))
            ax[i].set_ylim(0.0, ymax)
        ax[i].set_title(ttl)
        ax[i].grid(alpha=0.3)
        ax[i].set_xlabel("t (h)")

        ax2 = ax[i].twinx()
        ax2.step(tt[:-1], vv, where="post", lw=1.6, ls="--", color="tab:orange", alpha=0.85)
        ax2.set_ylabel(vlabel, color="tab:orange", fontsize=8)
        ax2.tick_params(axis="y", labelcolor="tab:orange", labelsize=8)

    for j in range(i + 1, len(ax)):
        ax[j].axis("off")
    plt.tight_layout()
    plt.show()

    flux_plot_series = [
        ("mu (1/h)", mu_hist),
        ("v_glu (mmol/gDW/h)", vglu_hist),
        ("v_fru (mmol/gDW/h)", vfru_hist),
        ("v_eth (mmol/gDW/h)", veth_hist),
        ("v_o2 (mmol/gDW/h)", vo2_hist),
        ("v_ATPM (mmol/gDW/h)", vatpm_hist),
        ("v_Prot (mmol/gDW/h)", vprot_hist),
        ("v_N (gN/gDW/h)", vn_eff_hist),
    ] + [
        (f"v_{aa_key} (mmol/gDW/h)", precursor_uptake_hist[aa_key]) for aa_key in aa_keys
    ] + [
        (f"v_turnover_{aa_key} (mmol/gDW/h)", aa_turnover_internal_hist[aa_key] / np.maximum(X_hist[:-1], EPS)) for aa_key in aa_keys
    ] + [
        (f"v_{akey} (mmol/gDW/h)", aroma_flux_hist[akey]) for akey in aroma_keys
    ]

    ncols_v, nrows_v = 4, math.ceil(len(flux_plot_series) / 4)
    fig_v, ax_v = plt.subplots(nrows_v, ncols_v, figsize=(5.2 * ncols_v, 3.2 * nrows_v), sharex=True, dpi=160)
    ax_v = np.atleast_1d(ax_v).ravel()
    for i, (ttl, vv) in enumerate(flux_plot_series):
        ax_v[i].step(tt[:-1], vv, where="post", lw=2, color="tab:orange")
        if np.isfinite(t_stationary):
            ax_v[i].axvline(t_stationary, ls="--", lw=1, color="tab:blue")
        ax_v[i].set_title(ttl)
        ax_v[i].grid(alpha=0.3)
        ax_v[i].set_xlabel("t (h)")
    for j in range(i + 1, len(ax_v)):
        ax_v[j].axis("off")
    plt.tight_layout()
    plt.show()

    return {
        "states_df": states_df,
        "flux_df": flux_df,
        "Y": Y,
        "mode_hist": mode_hist,
        "found_precursors": found_precursors,
        "found_aromas": found_aromas,
        "PROT_RXN_ID": PROT_RXN_ID,
        "v_interest_df": v_interest_df,
        "time_grid": tt,
        "dfba_dt": dt,
    }

res_turn = _run_dfba_turnover_aroma(
    scenario_name="Anaeróbico + turnover",
    aerobic_mode=False,
    o2_init_g_l=0.0,
)