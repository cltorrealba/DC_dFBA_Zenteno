import math
import re
import numpy as np
import pandas as pd
import cobra

# ======================================================================================
# dFBA SOA + pFBA con variable GAM (protein-turnover) + aromas
# ======================================================================================
import matplotlib.pyplot as plt

print("\n" + "=" * 100)
print("dFBA SOA + pFBA con GAM variable + protein-turnover + aromas")
print("=" * 100)

# ── Amino acid MW table (g/mol) ────────────────────────────────────────────────
AA_MW_TABLE = {
    's_0404': 89.09,   's_0542': 121.16,  's_0432': 133.11,  's_0748': 147.13,
    's_1314': 165.19,  's_0757': 75.07,   's_0832': 155.15,  's_0847': 131.17,
    's_1099': 146.19,  's_1077': 131.17,  's_1148': 149.21,  's_0430': 132.12,
    's_1379': 115.13,  's_0747': 146.14,  's_0428': 174.20,  's_1428': 105.09,
    's_1491': 119.12,  's_1561': 117.15,  's_1527': 204.23,  's_1533': 181.19,
}
CARB_MW_TABLE = {
    's_0001': 180.16,  's_0004': 180.16,  's_0509': 221.21,
    's_0773': 180.16,  's_1107': 180.16,  's_1520': 342.296,
}
RNA_MW_TABLE = {
    's_0423': 347.2212,  's_0526': 323.1965,
    's_0782': 363.22,    's_1545': 324.1813,
}

# ── Base GAM / NGAM ────────────────────────────────────────────────────────────
GAM_BASE   = 30.49    # mmol ATP/gDW  (base value in yeast-GEM biomass reaction)
NGAM_BASE  = 0.7      # mmol ATP/gDW/h  (already set via ATPM bounds)
RNA_FRAC   = 0.06     # g RNA / gDW  (fixed)

# ── protein content: dynamic variable (g/gDW) ─────────────────────────────────
PROT_CONTENT_0 = 0.46   # initial g protein / gDW  (Nissen 1997 yeast ≈ 0.46)
CARB_CONTENT_0 = 0.37   # initial g carbohydrate / gDW
K_DEATH        = 0.005
TURNOVER_LAMBDA = 0.03
VPROT_FLOOR_FRAC = 0.90
XA_FRACTION    = 1.0
N_AMMONIA_FRACTION = 0.50
N_TOTAL_DEPLETION_THRESHOLD = 1e-3
K_AA_UPTAKE_GROWTH = 0.08
BIOMASS_LOCK_FRAC  = 0.999
AA_UPTAKE_WEIGHT   = 1.0
EPS_LOC            = float(globals().get("EPS", 1e-9))

# ── cell biomass composition constants ──────────────────────────────────────────
GAM_COEFF_P = 16.965   # mmol ATP / (gProtein/gDW)
GAM_COEFF_R = 1.638    # mmol ATP / (gRNA/gDW)
GAM_COEFF_C = 5.210    # mmol ATP / (gCarb/gDW)

# ── helper: find biomass reaction and protein pseudo-rxn ───────────────────────
def _find_biomass_rxn(mdl):
    """Return biomass reaction (r_2111 or fallback)."""
    obj_id = globals().get("OBJ_ID", "r_2111")
    if obj_id in {r.id for r in mdl.reactions}:
        return mdl.reactions.get_by_id(obj_id)
    for r in mdl.reactions:
        if "biomass" in r.name.lower():
            return r
    raise RuntimeError("No se encontró reacción de biomasa.")

def _find_protein_pseudo_rxn(mdl):
    """Return the pseudo-rxn that produces 'protein' (r_4047 or search)."""
    pid = globals().get("PROT_RXN_ID", "r_4047")
    if pid in {r.id for r in mdl.reactions}:
        return mdl.reactions.get_by_id(pid)
    for r in mdl.reactions:
        if "protein" in r.name.lower():
            return r
    return None

def _strip_compartment(met_id):
    """Remove compartment suffix: 's_0404[c]' -> 's_0404'."""
    return re.sub(r'\[.*?\]$', '', met_id)

def _calculate_content(mdl):
    """
    Calculate Ptot (g/gDW), Ctot (g/gDW), Rtot (g/gDW) from current biomass stoichiometry.
    Mirrors Sanchez's calculateContent.
    """
    prot_rxn = _find_protein_pseudo_rxn(mdl)
    if prot_rxn is None:
        return PROT_CONTENT_0, CARB_CONTENT_0, RNA_FRAC

    # Find carbohydrate and RNA pseudo-rxns by metabolite name
    carb_rxn = None
    rna_rxn  = None
    for r in mdl.reactions:
        for m in r.metabolites:
            if m.name.lower() == "carbohydrate" and r.metabolites[m] == 1:
                carb_rxn = r
            if m.name.lower() == "rna" and r.metabolites[m] == 1:
                rna_rxn = r

    Ptot, Ctot, Rtot = 0.0, 0.0, 0.0
    for rxn, table, acc_name in [
        (prot_rxn, AA_MW_TABLE,   "Ptot"),
        (carb_rxn, CARB_MW_TABLE, "Ctot"),
        (rna_rxn,  RNA_MW_TABLE,  "Rtot"),
    ]:
        if rxn is None:
            continue
        for met, coef in rxn.metabolites.items():
            key = _strip_compartment(met.id)
            if key in table and coef != 1.0 and coef < 0:
                mw = table[key]
                val = abs(float(coef)) * mw / 1000.0   # mmol/gDW * g/mmol
                if acc_name == "Ptot":
                    Ptot += val
                elif acc_name == "Ctot":
                    Ctot += val
                else:
                    Rtot += val

    Ptot = Ptot if Ptot > 1e-6 else PROT_CONTENT_0
    Ctot = Ctot if Ctot > 1e-6 else CARB_CONTENT_0
    Rtot = Rtot if Rtot > 1e-6 else RNA_FRAC
    return Ptot, Ctot, Rtot


def _compute_full_gam(P, R, C, Pbase, Rbase, Cbase):
    """
    fullGAM = GAM + 16.965*Pfactor + 1.638*Rfactor + 5.210*Cfactor
    Mirrors Sanchez's changeGAM.
    """
    Pfactor = P / max(Pbase, EPS_LOC)
    Rfactor = R / max(Rbase, EPS_LOC)
    Cfactor = (Cbase + Pbase - P - R) / max(Cbase, EPS_LOC)
    Cfactor = max(0.0, Cfactor)
    full_gam = GAM_BASE + GAM_COEFF_P * Pfactor + GAM_COEFF_R * Rfactor + GAM_COEFF_C * Cfactor
    return float(full_gam)


def _apply_variable_gam(mtmp, full_gam):
    """
    Update stoichiometry of ATP/ADP/H2O/H+/phosphate in the biomass reaction
    to reflect full_gam (mmol/gDW).
    Mirrors the bioRxn update block in changeGAM.
    """
    bio_rxn = _find_biomass_rxn(mtmp)
    targets = {"atp", "adp", "h2o", "h+", "phosphate", "water", "proton"}
    updates = {}
    for met, coef in bio_rxn.metabolites.items():
        nm = met.name.lower()
        if any(t in nm for t in targets):
            updates[met] = float(np.sign(coef)) * round(full_gam, 4) - float(coef)
    if updates:
        bio_rxn.add_metabolites(updates)


# ── Get base composition once from the original model ───────────────────────────
Pbase_global, Cbase_global, Rbase_global = _calculate_content(model)
print(f"Base composition (gDW): Prot={Pbase_global:.4f}, Carb={Cbase_global:.4f}, RNA={Rbase_global:.4f}")
gam_base_check = _compute_full_gam(Pbase_global, Rbase_global, Cbase_global,
                                    Pbase_global, Rbase_global, Cbase_global)
print(f"GAM at base composition: {gam_base_check:.4f} mmol ATP/gDW")

# ── Re-use mappings from cell 19 if available ───────────────────────────────────
_found_precursors = globals().get("found_precursors", {})
_found_aromas     = globals().get("found_aromas", {})
_PROT_RXN_ID      = globals().get("PROT_RXN_ID", "r_4047")
_AA_MW            = globals().get("AA_MW", {"phe":0.16519,"leu":0.13117,"val":0.11715,"met":0.14921,"tyr":0.18119})
_AROMA_MW         = globals().get("AROMA_MW", {"pea":0.12217,"isoamyl":0.08815,"isobutanol":0.07412,"methionol":0.10619,"tyrosol":0.13816})
_N_atoms_map      = globals().get("N_atoms_map", {rid: 1.0 for rid in KINETIC_N_SOURCE_IDS})
_N_frac_map       = globals().get("N_frac_map", {rid: 0.02 for rid in KINETIC_N_SOURCE_IDS})

if not _found_precursors:
    raise RuntimeError("'found_precursors' not found. Run cell 19 first.")

aa_keys   = list(_found_precursors.keys())
aroma_keys = list(_found_aromas.keys())

AA_ALPHA_STATIC = float(globals().get("AA_ALPHA_STATIC", 1.0 / max(1, len(aa_keys))))
AA_ALPHA        = {k: AA_ALPHA_STATIC for k in aa_keys}

# Build tracked IDs
_required_ids = list(dict.fromkeys(
    [OBJ_ID, GLU_ID, FRU_ID, ETH_ID, O2_ID, ATPM_ID, _PROT_RXN_ID]
    + list(_found_precursors.values())
    + list(_found_aromas.values())
    + list(KINETIC_N_SOURCE_IDS)
))
_required_ids = [rid for rid in _required_ids if rid in rxn_by_id]

_base_lb = {rid: float(rxn_by_id[rid].lower_bound) for rid in _required_ids}
_base_ub = {rid: float(rxn_by_id[rid].upper_bound) for rid in _required_ids}

# ── Kinetic limits (reuse existing function if available) ───────────────────────
_kinetic_limits_fn = globals().get("_kinetic_limits")
if _kinetic_limits_fn is None:
    raise RuntimeError("'_kinetic_limits' not found. Run cell 15/17/19 first.")

N_VITAMIN_EX_IDS_LOC = globals().get("N_VITAMIN_EX_IDS", [])

# ── dFBA with variable GAM ────────────────────────────────────────────────────
def _run_dfba_variable_gam(
    *,
    scenario_name="Anaeróbico + GAM variable",
    aerobic_mode=False,
    o2_init_g_l=0.0,
    dfba_dt=1.0,
    t_end=72.0,
):
    dt   = float(dfba_dt)
    tt   = np.arange(0.0, t_end + 0.5 * dt, dt)
    n_st = len(tt) - 1

    N0_total   = float(globals().get("N0", 0.14))
    N0_ammonia = N0_total * N_AMMONIA_FRACTION
    N0_from_aa = max(0.0, N0_total - N0_ammonia)
    aa0_each   = (N0_from_aa / MW_N) / max(1, len(aa_keys))

    # State vector: [X, N_free, G, F, E, O2, Prot, Carb, AA_i..., Aroma_j...]
    n_aa    = len(aa_keys)
    n_aroma = len(aroma_keys)
    n_state = 8 + n_aa + n_aroma
    Y = np.zeros((len(tt), n_state), dtype=float)
    Y[0, 0] = 0.5
    Y[0, 1] = N0_ammonia
    Y[0, 2] = 110.0
    Y[0, 3] = 110.0
    Y[0, 4] = 0.0
    Y[0, 5] = o2_init_g_l
    Y[0, 6] = Y[0, 0] * PROT_CONTENT_0    # Prot [g/L]
    Y[0, 7] = Y[0, 0] * CARB_CONTENT_0    # Carb [g/L]
    for j in range(n_aa):
        Y[0, 8 + j] = aa0_each

    # Flux histories
    mu_h         = np.zeros(n_st)
    vglu_h       = np.zeros(n_st)
    vfru_h       = np.zeros(n_st)
    veth_h       = np.zeros(n_st)
    vo2_h        = np.zeros(n_st)
    vn_h         = np.zeros(n_st)
    vatpm_h      = np.zeros(n_st)
    vprot_h      = np.zeros(n_st)
    gam_h        = np.zeros(n_st)
    p_frac_h     = np.zeros(n_st)   # g protein / gDW
    c_frac_h     = np.zeros(n_st)   # g carb / gDW
    aroma_fh     = {k: np.zeros(n_st) for k in aroma_keys}
    aa_uptake_h  = {k: np.zeros(n_st) for k in aa_keys}
    aa_turn_h    = {k: np.zeros(n_st) for k in aa_keys}
    mode_h       = []

    for k in range(n_st):
        yk = np.maximum(Y[k, :], 0.0)
        cX, cN_free, cG, cF, cE, cO2, cProt, cCarb = yk[:8]
        aa_state = {ak: float(yk[8 + j]) for j, ak in enumerate(aa_keys)}

        aa_n_sum = MW_N * sum(aa_state.values())
        cN_total = cN_free + aa_n_sum

        # Dynamic protein / carb fractions
        P_frac = cProt / max(cX, EPS_LOC)   # g/gDW
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

            # Reset bounds
            for rid in _required_ids:
                rb[rid].lower_bound = _base_lb[rid]
                rb[rid].upper_bound = _base_ub[rid]

            # Sugars
            if GLU_ID in rb:
                rb[GLU_ID].lower_bound = max(rb[GLU_ID].lower_bound, -lim_glu)
                rb[GLU_ID].upper_bound = min(rb[GLU_ID].upper_bound, 0.0)
            if FRU_ID in rb:
                rb[FRU_ID].lower_bound = max(rb[FRU_ID].lower_bound, -lim_fru)
                rb[FRU_ID].upper_bound = min(rb[FRU_ID].upper_bound, 0.0)

            # N kinetic
            for rid in KINETIC_N_SOURCE_IDS:
                if rid in rb:
                    Li = float(lim_n_map.get(rid, 0.0))
                    rb[rid].lower_bound = max(rb[rid].lower_bound, -Li)
                    rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)

            # Vitamins/N depletion
            if cN_total <= N_TOTAL_DEPLETION_THRESHOLD:
                for rid in N_VITAMIN_EX_IDS_LOC:
                    if rid in rb:
                        rb[rid].lower_bound = max(float(rb[rid].lower_bound), 0.0)
                        rb[rid].upper_bound = min(float(rb[rid].upper_bound), 0.0)

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

            # Apply variable GAM to biomass stoichiometry
            _apply_variable_gam(mtmp, full_gam)

            # AA uptake caps
            turnover_caps = {}
            if cN_total > N_TOTAL_DEPLETION_THRESHOLD:
                for ak, rid in _found_precursors.items():
                    if rid in rb:
                        q_i = K_AA_UPTAKE_GROWTH * aa_state[ak] / max(cX, EPS_LOC)
                        q_i = max(0.0, q_i)
                        rb[rid].lower_bound = -q_i
                        rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)
                        turnover_caps[ak] = q_i
            else:
                for rid in KINETIC_N_SOURCE_IDS:
                    if rid in rb:
                        rb[rid].lower_bound = max(float(rb[rid].lower_bound), 0.0)
                        rb[rid].upper_bound = min(float(rb[rid].upper_bound), 0.0)
                Xa = XA_FRACTION * cX
                for ak, rid in _found_precursors.items():
                    if rid in rb:
                        q_i = TURNOVER_LAMBDA * cProt * AA_ALPHA[ak] / max(Xa, EPS_LOC)
                        q_i = max(0.0, q_i)
                        rb[rid].lower_bound = -q_i
                        rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)
                        turnover_caps[ak] = q_i

            # Phase switch
            if cN_total > N_TOTAL_DEPLETION_THRESHOLD:
                mode_h.append("BIOMASS")
                mtmp.objective = OBJ_ID
                mtmp.objective_direction = "max"
                sol_mu = mtmp.optimize()
                if getattr(sol_mu, "status", "") != "optimal":
                    mode_h[-1] = "NO_OPT"; continue
                mu_star = max(0.0, float(sol_mu.fluxes.get(OBJ_ID, 0.0)))
                rb[OBJ_ID].lower_bound = max(float(rb[OBJ_ID].lower_bound),
                                              BIOMASS_LOCK_FRAC * mu_star)

                # Secondary obj: maximise AA uptake
                uptake_obj = {mtmp.reactions.get_by_id(rid): -AA_UPTAKE_WEIGHT
                              for rid in _found_precursors.values() if rid in rb}
                mtmp.objective = uptake_obj
                mtmp.objective_direction = "max"
                sol = cobra.flux_analysis.pfba(mtmp, fraction_of_optimum=1.0)
                if getattr(sol, "status", "") != "optimal":
                    mode_h[-1] = "NO_OPT"; continue
            else:
                mode_h.append("TURNOVER_ATPM")
                rb[OBJ_ID].lower_bound = 0.0
                rb[OBJ_ID].upper_bound = 0.0
                rb[ATPM_ID].lower_bound = max(float(rb[ATPM_ID].lower_bound),
                                               ATPM_LB_NO_GROWTH)
                rb[ATPM_ID].upper_bound = ATPM_UB_NO_GROWTH
                mtmp.objective = ATPM_ID
                mtmp.objective_direction = "max"
                sol = cobra.flux_analysis.pfba(mtmp, fraction_of_optimum=1.0)
                if getattr(sol, "status", "") != "optimal":
                    mode_h[-1] = "NO_OPT"; continue

            v_obj  = float(sol.fluxes.get(OBJ_ID,   0.0))
            v_glu  = float(sol.fluxes.get(GLU_ID,   0.0))
            v_fru  = float(sol.fluxes.get(FRU_ID,   0.0))
            v_eth  = float(sol.fluxes.get(ETH_ID,   0.0))
            v_o2   = float(sol.fluxes.get(O2_ID,    0.0))
            v_atpm = float(sol.fluxes.get(ATPM_ID,  0.0))
            v_prot = float(sol.fluxes.get(_PROT_RXN_ID, 0.0))

            if cN_total > N_TOTAL_DEPLETION_THRESHOLD:
                vn_eff = 0.0
                for rid in KINETIC_N_SOURCE_IDS:
                    if rid in _base_lb:
                        vi_ = float(sol.fluxes.get(rid, 0.0))
                        vn_eff += max(0.0, -vi_) * _N_atoms_map.get(rid, 1.0) * MW_N
            else:
                vn_eff = 0.0

            for j, ak in enumerate(aa_keys):
                rid = _found_precursors[ak]
                v_aa = float(sol.fluxes.get(rid, 0.0))
                aa_uptake_h[ak][k] = max(0.0, -v_aa)
                aa_turn_h[ak][k]   = (TURNOVER_LAMBDA * cProt * AA_ALPHA[ak]
                                       if mode_h[-1] == "TURNOVER_ATPM" else 0.0)

            for ak in aroma_keys:
                rid = _found_aromas[ak]
                aroma_fh[ak][k] = max(0.0, float(sol.fluxes.get(rid, 0.0)))

        mu_h[k]    = max(0.0, v_obj)
        vglu_h[k]  = max(0.0, -v_glu)
        vfru_h[k]  = max(0.0, -v_fru)
        veth_h[k]  = max(0.0, v_eth)
        vo2_h[k]   = max(0.0, -v_o2)
        vn_h[k]    = max(0.0, vn_eff)
        vatpm_h[k] = max(0.0, v_atpm)
        vprot_h[k] = max(0.0, v_prot)

        # ODE integration (Euler)
        dX    = v_obj * cX if mode_h[-1] == "BIOMASS" else 0.0
        dN    = -vn_eff * cX if mode_h[-1] == "BIOMASS" else 0.0
        dG    = -MW_GLU * max(0.0, -v_glu) * cX
        dF    = -MW_FRU * max(0.0, -v_fru) * cX
        dE    =  MW_ETH * max(0.0,  v_eth) * cX
        dO2   = -MW_O2  * max(0.0, -v_o2)  * cX

        Xa  = XA_FRACTION * cX
        dProt = (
            max(0.0, v_obj) * cX * PROT_CONTENT_0
            - cProt * K_DEATH
            + Xa * max(0.0, v_prot)
            - TURNOVER_LAMBDA * cProt
        )
        # Carb: simple decrease proportional to biomass growth (synthesis - dilution)
        dCarb = max(0.0, v_obj) * cX * CARB_CONTENT_0 - cCarb * K_DEATH

        Y[k+1, 0] = max(0.0, cX    + dt * dX)
        Y[k+1, 1] = max(0.0, cN_free + dt * dN)
        Y[k+1, 2] = max(0.0, cG    + dt * dG)
        Y[k+1, 3] = max(0.0, cF    + dt * dF)
        Y[k+1, 4] = max(0.0, cE    + dt * dE)
        Y[k+1, 5] = max(0.0, cO2   + dt * dO2)
        Y[k+1, 6] = max(0.0, cProt + dt * dProt)
        Y[k+1, 7] = max(0.0, cCarb + dt * dCarb)

        for j, ak in enumerate(aa_keys):
            v_aa_i = float(sol.fluxes.get(_found_precursors[ak], 0.0))
            dAA = cX * v_aa_i if mode_h[-1] == "BIOMASS" else 0.0
            Y[k+1, 8+j] = max(0.0, Y[k, 8+j] + dt * dAA)

        for j, ak in enumerate(aroma_keys):
            mw = _AROMA_MW.get(ak, 0.1)
            dA = mw * aroma_fh[ak][k] * cX
            Y[k+1, 8+n_aa+j] = max(0.0, Y[k, 8+n_aa+j] + dt * dA)

    # Unpack states
    X_h    = Y[:, 0]
    Nf_h   = Y[:, 1]
    G_h    = Y[:, 2]
    F_h    = Y[:, 3]
    E_h    = Y[:, 4]
    O2_h   = Y[:, 5]
    Prot_h = Y[:, 6]
    Carb_h = Y[:, 7]
    aa_conc_h  = {ak: Y[:, 8+j]         for j, ak in enumerate(aa_keys)}
    aroma_conc_h = {ak: Y[:, 8+n_aa+j]  for j, ak in enumerate(aroma_keys)}
    aroma_mg_h   = {ak: 1000.0 * aroma_conc_h[ak] for ak in aroma_keys}

    aa_n_h = MW_N * np.sum(
        np.column_stack([aa_conc_h[ak] for ak in aa_keys]), axis=1)
    N_h    = Nf_h + aa_n_h

    # Protein/carb fraction vs X
    P_frac_h_full = Prot_h / np.maximum(X_h, EPS_LOC)
    C_frac_h_full = Carb_h / np.maximum(X_h, EPS_LOC)

    # GAM at every time point (for plot)
    gam_full = np.array([
        _compute_full_gam(float(P_frac_h_full[i]), RNA_FRAC, float(C_frac_h_full[i]),
                          Pbase_global, Rbase_global, Cbase_global)
        for i in range(len(tt))
    ])

    t_stat = float(tt[np.argmax(np.array([m == "TURNOVER_ATPM" for m in mode_h]))]) \
        if any(m == "TURNOVER_ATPM" for m in mode_h) else np.nan
    biomass_steps   = sum(m == "BIOMASS"       for m in mode_h)
    turnover_steps  = sum(m == "TURNOVER_ATPM" for m in mode_h)

    print(f"\n[{scenario_name}]")
    print(f"  Pasos BIOMASS={biomass_steps}, TURNOVER_ATPM={turnover_steps}")
    print(f"  t_stationary ≈ {t_stat:.1f} h")
    print(f"  X_final={X_h[-1]:.4f} gDW/L  |  E_final={E_h[-1]:.4f} g/L  |  Prot_final/X={P_frac_h_full[-1]:.4f}")
    print(f"  GAM_initial={gam_full[0]:.3f}  →  GAM_final={gam_full[-1]:.3f}  mmol ATP/gDW")

    # ── States DataFrame ──────────────────────────────────────────────────────
    states_df = pd.DataFrame({
        "t_h": tt, "X_gDW_L": X_h, "N_gN_L": N_h, "N_free_gN_L": Nf_h,
        "G_g_L": G_h, "F_g_L": F_h, "E_g_L": E_h, "O2_g_L": O2_h,
        "Prot_g_L": Prot_h, "Carb_g_L": Carb_h,
        "P_frac_gProt_gDW": P_frac_h_full, "C_frac_gCarb_gDW": C_frac_h_full,
        "GAM_mmol_gDW": gam_full,
    })
    for ak in aa_keys:
        states_df[f"AA_{ak}_mmol_L"] = aa_conc_h[ak]
    for ak in aroma_keys:
        states_df[f"{ak}_g_L"]  = aroma_conc_h[ak]
        states_df[f"{ak}_mg_L"] = aroma_mg_h[ak]

    # ── Flux DataFrame ────────────────────────────────────────────────────────
    flux_df = pd.DataFrame({
        "t_h": tt[:-1], "mode": mode_h,
        "mu_h_inv": mu_h, "v_glu": vglu_h, "v_fru": vfru_h,
        "v_eth": veth_h, "v_o2": vo2_h, "v_N_uptake": vn_h,
        "v_ATPM": vatpm_h, "v_Prot": vprot_h, "GAM": gam_h,
    })
    for ak in aa_keys:
        flux_df[f"v_uptake_{ak}"] = aa_uptake_h[ak]
        flux_df[f"v_turn_{ak}"]   = aa_turn_h[ak]
    for ak in aroma_keys:
        flux_df[f"v_{ak}"] = aroma_fh[ak]

    display(flux_df)

    # ═══════════════════════════════════════════════════════════════════════════
    # PLOTS
    # ═══════════════════════════════════════════════════════════════════════════

    # ── 1. Main states ────────────────────────────────────────────────────────
    state_series = [
        ("X (gDW/L)",   X_h,    mu_h,   "mu (1/h)"),
        ("N (gN/L)",    N_h,    vn_h,   "v_N (gN/gDW/h)"),
        ("Glucosa (g/L)", G_h,  vglu_h, "v_glu"),
        ("Fructosa (g/L)", F_h, vfru_h, "v_fru"),
        ("Etanol (g/L)", E_h,   veth_h, "v_eth"),
        ("O2 (g/L)",    O2_h,   vo2_h,  "v_o2"),
        ("Prot (g/L)",  Prot_h, vprot_h,"v_Prot"),
        ("Carb (g/L)",  Carb_h, np.zeros(n_st), "—"),
    ] + [
        (f"AA_{ak} (mmol/L)", aa_conc_h[ak], aa_uptake_h[ak], f"v_{ak}")
        for ak in aa_keys
    ] + [
        (f"{ak} (mg/L)", aroma_mg_h[ak], aroma_fh[ak], f"v_{ak}")
        for ak in aroma_keys
    ]

    ncols, nrows = 4, math.ceil(len(state_series) / 4)
    fig, ax = plt.subplots(nrows, ncols, figsize=(5.2*ncols, 3.5*nrows), sharex=True, dpi=120)
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

    # ── 2. GAM evolution ──────────────────────────────────────────────────────
    fig_gam, ax_gam = plt.subplots(1, 3, figsize=(15, 4), dpi=120)

    ax_gam[0].plot(tt, gam_full, lw=2.5, color="tab:red")
    ax_gam[0].axhline(GAM_BASE, ls="--", lw=1.2, color="gray", label=f"GAM_base={GAM_BASE}")
    if np.isfinite(t_stat):
        ax_gam[0].axvline(t_stat, ls=":", lw=1.5, color="navy", alpha=0.6, label=f"t_stat={t_stat:.0f}h")
    ax_gam[0].set_title("GAM (mmol ATP/gDW)", fontsize=10)
    ax_gam[0].set_xlabel("t (h)")
    ax_gam[0].legend(fontsize=8)
    ax_gam[0].grid(alpha=0.3)

    ax_gam[1].plot(tt, P_frac_h_full, lw=2.5, color="tab:green", label="Protein (g/gDW)")
    ax_gam[1].plot(tt, C_frac_h_full, lw=2.5, color="tab:purple", label="Carb (g/gDW)")
    ax_gam[1].axhline(Pbase_global, ls="--", lw=1, color="tab:green", alpha=0.5)
    ax_gam[1].axhline(Cbase_global, ls="--", lw=1, color="tab:purple", alpha=0.5)
    if np.isfinite(t_stat):
        ax_gam[1].axvline(t_stat, ls=":", lw=1.5, color="navy", alpha=0.6)
    ax_gam[1].set_title("Composición celular (g/gDW)", fontsize=10)
    ax_gam[1].set_xlabel("t (h)")
    ax_gam[1].legend(fontsize=8)
    ax_gam[1].grid(alpha=0.3)

    ax_gam[2].plot(tt, Prot_h, lw=2.5, color="tab:green",  label="Prot (g/L)")
    ax_gam[2].plot(tt, X_h,    lw=2.5, color="tab:blue",   label="X (gDW/L)")
    ax_gam[2].plot(tt, Carb_h, lw=2.5, color="tab:purple", label="Carb (g/L)")
    if np.isfinite(t_stat):
        ax_gam[2].axvline(t_stat, ls=":", lw=1.5, color="navy", alpha=0.6)
    ax_gam[2].set_title("Biomasa, Proteína y Carbohidrato (g/L)", fontsize=10)
    ax_gam[2].set_xlabel("t (h)")
    ax_gam[2].legend(fontsize=8)
    ax_gam[2].grid(alpha=0.3)

    plt.suptitle(f"GAM variable y contenido proteico – {scenario_name}", fontsize=11)
    plt.tight_layout()
    plt.show()

    # ── 3. Aroma profiles ────────────────────────────────────────────────────
    if aroma_keys:
        ncols_a = min(len(aroma_keys), 5)
        nrows_a = math.ceil(len(aroma_keys) / ncols_a)
        fig_a, ax_a = plt.subplots(nrows_a, ncols_a,
                                    figsize=(4.5*ncols_a, 3.5*nrows_a), dpi=120, squeeze=False)
        ax_a = ax_a.ravel()
        for i, ak in enumerate(aroma_keys):
            ax_a[i].plot(tt, aroma_mg_h[ak], lw=2.5, color="tab:brown")
            ax_a[i].set_title(f"{ak} (mg/L)", fontsize=9)
            ax_a[i].set_xlabel("t (h)")
            ax_a[i].grid(alpha=0.3)
            if np.isfinite(t_stat):
                ax_a[i].axvline(t_stat, ls=":", lw=1.2, color="navy", alpha=0.5)
        for j in range(i+1, len(ax_a)):
            ax_a[j].axis("off")
        plt.suptitle(f"Aromas fermentativos – {scenario_name}", fontsize=11)
        plt.tight_layout()
        plt.show()

    # Save CSV
    sfx = "aerobic_vargam" if aerobic_mode else "anaerobic_vargam"
    out_file = OUT_DIR / f"dfba_vargam_states_{sfx}.csv"
    states_df.to_csv(out_file, index=False)
    print(f"Archivo guardado: {out_file}")

    return {
        "states_df":    states_df,
        "flux_df":      flux_df,
        "Y":            Y,
        "mode_h":       mode_h,
        "gam_full":     gam_full,
        "P_frac_full":  P_frac_h_full,
        "C_frac_full":  C_frac_h_full,
        "time_grid_vg": tt,
    }


# ── Run anaerobic scenario ─────────────────────────────────────────────────────
res_vargam_ana = _run_dfba_variable_gam(
    scenario_name="Anaeróbico + GAM variable",
    aerobic_mode=True,
    o2_init_g_l=0.008,
    dfba_dt=1.0,
    t_end=72.0,
)