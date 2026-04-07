"""
=============================================================================
EXPORT SCRIPT: Generate dFVB inputs from Notebook 01 model state
=============================================================================
Run this script AS A CELL in Notebook 01 (after executing all cells up to 46)
so that the model object and all runtime variables are in scope.

Outputs generated in Material_dFVB/:
  1. S_matrix.csv            - Stoichiometric matrix (m x n)
  2. reaction_info.csv       - Reaction metadata (id, name, bounds, subsystem)
  3. metabolite_info.csv     - Metabolite metadata (id, name, compartment, formula)
  4. bounds_static.csv       - Base static bounds (lb, ub) for each reaction
  5. bounds_dynamic_rules.csv- Rules for dynamic bounds (Zenteno kinetics)
  6. objective_config.csv    - Objective functions and switching criteria
  7. reaction_indices_of_interest.csv - Indices of aroma, AA, key reactions
  8. initial_conditions.csv  - y0 vector
  9. kinetic_parameters.csv  - Zenteno ODE kinetic parameters
  10. ode_state_definitions.csv - State vector structure
  11. gam_parameters.csv     - Variable GAM parameters
  12. nitrogen_config.csv    - N source mapping and atoms
  13. aroma_config.csv       - Aroma exchange IDs and MW
  14. aa_precursor_config.csv- AA precursor exchange IDs and MW
  15. turnover_parameters.csv- Protein turnover / stationary phase parameters
  16. tracked_bounds_log.csv - All set_bounds_logged entries
=============================================================================
"""

import os
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import sparse
import cobra

# ────────────────────────────────────────────────────────────────────────────
# Destination folder
# ────────────────────────────────────────────────────────────────────────────
OUT = Path("Material_dFVB") / "csv_exports"
OUT.mkdir(parents=True, exist_ok=True)
print(f"[dFVB Export] Output → {OUT.resolve()}")

# ────────────────────────────────────────────────────────────────────────────
# 1. Stoichiometric matrix  (S_matrix.csv)
# ────────────────────────────────────────────────────────────────────────────
print("[1/16] Stoichiometric matrix ...")
S = cobra.util.array.create_stoichiometric_matrix(model, array_type="dense")
rxn_ids = [r.id for r in model.reactions]
met_ids = [m.id for m in model.metabolites]

S_df = pd.DataFrame(S, index=met_ids, columns=rxn_ids)
S_df.index.name = "metabolite_id"
S_df.to_csv(OUT / "S_matrix.csv")
print(f"  S_matrix.csv: {S.shape[0]} metabolites x {S.shape[1]} reactions")

# Also save sparse triplet format (more MATLAB-friendly for large S)
rows, cols = np.nonzero(S)
vals = S[rows, cols]
sparse_df = pd.DataFrame({
    "row": rows,
    "col": cols,
    "value": vals,
    "metabolite_id": [met_ids[r] for r in rows],
    "reaction_id": [rxn_ids[c] for c in cols],
})
sparse_df.to_csv(OUT / "S_matrix_sparse.csv", index=False)
print(f"  S_matrix_sparse.csv: {len(vals)} non-zero entries")

# ────────────────────────────────────────────────────────────────────────────
# 2. Reaction info  (reaction_info.csv)
# ────────────────────────────────────────────────────────────────────────────
print("[2/16] Reaction info ...")
rxn_rows = []
for i, r in enumerate(model.reactions):
    rxn_rows.append({
        "index_1based": i + 1,
        "reaction_id": r.id,
        "reaction_name": r.name,
        "lower_bound": r.lower_bound,
        "upper_bound": r.upper_bound,
        "subsystem": r.subsystem,
        "is_exchange": r.id in {ex.id for ex in model.exchanges},
        "gene_reaction_rule": r.gene_reaction_rule,
        "reaction_string": r.reaction,
    })
pd.DataFrame(rxn_rows).to_csv(OUT / "reaction_info.csv", index=False)
print(f"  reaction_info.csv: {len(rxn_rows)} reactions")

# ────────────────────────────────────────────────────────────────────────────
# 3. Metabolite info  (metabolite_info.csv)
# ────────────────────────────────────────────────────────────────────────────
print("[3/16] Metabolite info ...")
met_rows = []
for i, m in enumerate(model.metabolites):
    met_rows.append({
        "index_1based": i + 1,
        "metabolite_id": m.id,
        "metabolite_name": m.name,
        "compartment": m.compartment,
        "formula": getattr(m, "formula", ""),
        "charge": getattr(m, "charge", ""),
    })
pd.DataFrame(met_rows).to_csv(OUT / "metabolite_info.csv", index=False)
print(f"  metabolite_info.csv: {len(met_rows)} metabolites")

# ────────────────────────────────────────────────────────────────────────────
# 4. Static bounds  (bounds_static.csv)
# ────────────────────────────────────────────────────────────────────────────
print("[4/16] Static bounds ...")
bounds_rows = []
for i, r in enumerate(model.reactions):
    bounds_rows.append({
        "index_1based": i + 1,
        "reaction_id": r.id,
        "lb": r.lower_bound,
        "ub": r.upper_bound,
    })
pd.DataFrame(bounds_rows).to_csv(OUT / "bounds_static.csv", index=False)

# Also export as pure numeric vectors (MATLAB-friendly)
lb_vec = np.array([r.lower_bound for r in model.reactions])
ub_vec = np.array([r.upper_bound for r in model.reactions])
np.savetxt(OUT / "lb_vector.csv", lb_vec, delimiter=",", fmt="%.10g")
np.savetxt(OUT / "ub_vector.csv", ub_vec, delimiter=",", fmt="%.10g")
print(f"  bounds_static.csv, lb_vector.csv, ub_vector.csv")

# ────────────────────────────────────────────────────────────────────────────
# 5. Dynamic bounds rules  (bounds_dynamic_rules.csv)
# ────────────────────────────────────────────────────────────────────────────
print("[5/16] Dynamic bounds rules ...")
_EPS = float(globals().get("EPS", 1e-9))
_O2_VMAX = float(globals().get("O2_VMAX_UPTAKE", 2.0))
_O2_KO = float(globals().get("O2_KO_G_L", 0.002))
_O2_DEPL = float(globals().get("O2_DEPLETION_THRESHOLD", 1e-6))
_APPLY_CAPS = bool(globals().get("APPLY_PRODUCT_CAPS", False))

dynamic_rules = [
    {"reaction_id": "r_1714", "variable": "glucose",   "rule": "lb = max(base_lb, -lim_glu); ub = min(base_ub, 0)", "kinetics": "Monod(cG, Kg0) * EtOH_inhib(cE, Kie0)"},
    {"reaction_id": "r_1709", "variable": "fructose",  "rule": "lb = max(base_lb, -lim_fru); ub = min(base_ub, 0)", "kinetics": "Monod(cF, Kf0) * Glc_inhib(cG, Kig0) * EtOH_inhib(cE, Kie0)"},
    {"reaction_id": "r_1992", "variable": "O2",        "rule": "anaerobic: lb=0, ub=0; aerobic: lb=-lim_o2, ub=0", "kinetics": f"Monod(cO2, {_O2_KO}) * {_O2_VMAX}"},
    {"reaction_id": "r_1761", "variable": "ethanol",   "rule": f"APPLY_PRODUCT_CAPS={_APPLY_CAPS}: ub=min(base_ub, lim_eth)", "kinetics": "from sugar yields"},
    {"reaction_id": "r_2111", "variable": "biomass",   "rule": f"BIOMASS mode: lb=LOCK_FRAC*mu*; TURNOVER: lb=ub=0", "kinetics": "Monod(cN, Kn0)"},
    {"reaction_id": "r_4046", "variable": "ATPM",      "rule": "BIOMASS: lb=0.7, ub=1000; TURNOVER: lb=0.7, ub=1000", "kinetics": "fixed bounds"},
]
# Add N sources
for rid in KINETIC_N_SOURCE_IDS:
    dynamic_rules.append({
        "reaction_id": rid,
        "variable": f"N_source ({rxn_by_id[rid].name if rid in rxn_by_id else rid})",
        "rule": "BIOMASS: lb=max(base_lb, -Li), ub=min(base_ub,0); TURNOVER: lb=ub=0",
        "kinetics": "Li = vN_total * N_frac_map[rid]",
    })
# Add AA precursors
for ak, rid in (globals().get("found_precursors", {}) or {}).items():
    dynamic_rules.append({
        "reaction_id": rid,
        "variable": f"AA_precursor_{ak} ({rxn_by_id[rid].name if rid in rxn_by_id else rid})",
        "rule": "BIOMASS: lb=-K_AA*aa_i/X, ub=0; TURNOVER: lb=-lambda*Prot*alpha/Xa, ub=0",
        "kinetics": f"K_AA_UPTAKE_GROWTH={K_AA_UPTAKE_GROWTH}, TURNOVER_LAMBDA={TURNOVER_LAMBDA}",
    })

pd.DataFrame(dynamic_rules).to_csv(OUT / "bounds_dynamic_rules.csv", index=False)

# ────────────────────────────────────────────────────────────────────────────
# 6. Objective config  (objective_config.csv)
# ────────────────────────────────────────────────────────────────────────────
print("[6/16] Objective configuration ...")
_LOCK = float(globals().get("BIOMASS_LOCK_FRAC", 0.999))
_VPROT_FLOOR = float(globals().get("VPROT_FLOOR_FRAC", 0.90))
_AA_W = float(globals().get("AA_UPTAKE_WEIGHT", 1.0))

obj_rows = [
    {
        "phase": "BIOMASS",
        "step": "A",
        "description": "Maximize growth rate (mu)",
        "objective_reaction": "r_2111",
        "direction": "max",
        "constraint": "none",
    },
    {
        "phase": "BIOMASS",
        "step": "B",
        "description": f"Lock growth: mu >= {_LOCK} * mu_star",
        "objective_reaction": "r_2111",
        "direction": "constraint",
        "constraint": f"lb(r_2111) = {_LOCK} * mu_star",
    },
    {
        "phase": "BIOMASS",
        "step": "C",
        "description": f"Secondary: maximize AA uptake (w={_AA_W}) + pFBA",
        "objective_reaction": "AA_precursors (negative sum)",
        "direction": "max",
        "constraint": f"pFBA at fraction_of_optimum=1.0",
    },
    {
        "phase": "TURNOVER_ATPM",
        "step": "A",
        "description": "Close growth: r_2111 lb=ub=0",
        "objective_reaction": "r_2111",
        "direction": "fix",
        "constraint": "lb=0, ub=0",
    },
    {
        "phase": "TURNOVER_ATPM",
        "step": "B",
        "description": f"Maximize protein production, lock at {_VPROT_FLOOR}*vprot_star",
        "objective_reaction": "r_4047",
        "direction": "max then lock",
        "constraint": f"lb(r_4047) = {_VPROT_FLOOR} * vprot_star",
    },
    {
        "phase": "TURNOVER_ATPM",
        "step": "C",
        "description": "Maximize ATPM + pFBA",
        "objective_reaction": "r_4046",
        "direction": "max",
        "constraint": "pFBA at fraction_of_optimum=1.0",
    },
]
obj_df = pd.DataFrame(obj_rows)
obj_df.to_csv(OUT / "objective_config.csv", index=False)

# Phase switching
switch_df = pd.DataFrame([{
    "criterion": "N_total > N_TOTAL_DEPLETION_THRESHOLD",
    "N_TOTAL_DEPLETION_THRESHOLD_gN_L": float(globals().get("N_TOTAL_DEPLETION_THRESHOLD", 1e-3)),
    "formula": "N_total = N_free + MW_N * sum(AA_i)",
    "MW_N_g_per_mmol": float(globals().get("MW_N", 0.014007)),
    "if_true": "BIOMASS",
    "if_false": "TURNOVER_ATPM",
}])
switch_df.to_csv(OUT / "phase_switching_criterion.csv", index=False)

# ────────────────────────────────────────────────────────────────────────────
# 7. Reaction indices of interest  (reaction_indices_of_interest.csv)
# ────────────────────────────────────────────────────────────────────────────
print("[7/16] Reaction indices of interest ...")
rxn_id_to_idx = {r.id: i + 1 for i, r in enumerate(model.reactions)}

interest_rows = []
# Core reactions
core_map = {
    "OBJ_ID (biomass)": globals().get("OBJ_ID", "r_2111"),
    "GLU_ID (glucose)": globals().get("GLU_ID", "r_1714"),
    "FRU_ID (fructose)": globals().get("FRU_ID", "r_1709"),
    "ETH_ID (ethanol)": globals().get("ETH_ID", "r_1761"),
    "O2_ID (oxygen)": globals().get("O2_ID", "r_1992"),
    "ATPM_ID (maintenance)": globals().get("ATPM_ID", "r_4046"),
    "PROT_RXN_ID (protein)": globals().get("PROT_RXN_ID", "r_4047"),
}
for label, rid in core_map.items():
    idx = rxn_id_to_idx.get(rid, -1)
    rname = model.reactions.get_by_id(rid).name if rid in rxn_id_to_idx else ""
    interest_rows.append({
        "category": "core",
        "label": label,
        "reaction_id": rid,
        "reaction_name": rname,
        "index_1based": idx,
    })

# Aroma reactions
for ak, rid in (globals().get("found_aromas", {}) or {}).items():
    idx = rxn_id_to_idx.get(rid, -1)
    rname = model.reactions.get_by_id(rid).name if rid in rxn_id_to_idx else ""
    mw = (globals().get("AROMA_MW", {}) or {}).get(ak, "")
    interest_rows.append({
        "category": "aroma",
        "label": f"aroma_{ak}",
        "reaction_id": rid,
        "reaction_name": rname,
        "index_1based": idx,
        "MW_g_per_mmol": mw,
    })

# AA precursor reactions
for ak, rid in (globals().get("found_precursors", {}) or {}).items():
    idx = rxn_id_to_idx.get(rid, -1)
    rname = model.reactions.get_by_id(rid).name if rid in rxn_id_to_idx else ""
    mw = (globals().get("AA_MW", {}) or {}).get(ak, "")
    interest_rows.append({
        "category": "aa_precursor",
        "label": f"AA_{ak}",
        "reaction_id": rid,
        "reaction_name": rname,
        "index_1based": idx,
        "MW_g_per_mmol": mw,
    })

# Kinetic N sources
for rid in (globals().get("KINETIC_N_SOURCE_IDS", []) or []):
    idx = rxn_id_to_idx.get(rid, -1)
    rname = model.reactions.get_by_id(rid).name if rid in rxn_id_to_idx else ""
    n_atoms = (globals().get("N_atoms_map", {}) or {}).get(rid, 1.0)
    interest_rows.append({
        "category": "nitrogen_source",
        "label": f"N_source_{rid}",
        "reaction_id": rid,
        "reaction_name": rname,
        "index_1based": idx,
        "N_atoms": n_atoms,
    })

pd.DataFrame(interest_rows).to_csv(OUT / "reaction_indices_of_interest.csv", index=False)
print(f"  {len(interest_rows)} reactions of interest")

# ────────────────────────────────────────────────────────────────────────────
# 8. Initial conditions  (initial_conditions.csv)
# ────────────────────────────────────────────────────────────────────────────
print("[8/16] Initial conditions ...")
_MW_N = float(globals().get("MW_N", 0.014007))
_N0 = float(globals().get("N0", 0.14))
_N_AMM_FRAC = float(globals().get("N_AMMONIA_FRACTION", 0.5))
_N0_ammonia = _N0 * _N_AMM_FRAC
_N0_from_aa = max(0.0, _N0 - _N0_ammonia)
_aa_keys = list((globals().get("found_precursors", {}) or {}).keys())
_aroma_keys = list((globals().get("found_aromas", {}) or {}).keys())
_n_aa = len(_aa_keys)
_n_aroma = len(_aroma_keys)
_aa0_each = (_N0_from_aa / _MW_N) / max(1, _n_aa) if _n_aa > 0 else 0.0

ic_rows = [
    {"state_index": 1, "state_name": "X", "description": "Biomass", "unit": "gDW/L", "initial_value": 0.5},
    {"state_index": 2, "state_name": "N_free", "description": "Free nitrogen (NH4+)", "unit": "gN/L", "initial_value": _N0_ammonia},
    {"state_index": 3, "state_name": "G", "description": "Glucose", "unit": "g/L", "initial_value": 110.0},
    {"state_index": 4, "state_name": "F", "description": "Fructose", "unit": "g/L", "initial_value": 110.0},
    {"state_index": 5, "state_name": "E", "description": "Ethanol", "unit": "g/L", "initial_value": 0.0},
    {"state_index": 6, "state_name": "O2", "description": "Dissolved oxygen", "unit": "g/L", "initial_value": 0.0},
    {"state_index": 7, "state_name": "Prot", "description": "Intracellular protein", "unit": "g/L", "initial_value": 0.5 * 0.46},
    {"state_index": 8, "state_name": "Carb", "description": "Intracellular carbohydrate", "unit": "g/L", "initial_value": 0.5 * 0.37},
]
for j, ak in enumerate(_aa_keys):
    ic_rows.append({
        "state_index": 9 + j,
        "state_name": f"AA_{ak}",
        "description": f"Amino acid {ak} in medium",
        "unit": "mmol/L",
        "initial_value": _aa0_each,
    })
for j, ak in enumerate(_aroma_keys):
    ic_rows.append({
        "state_index": 9 + _n_aa + j,
        "state_name": f"Aroma_{ak}",
        "description": f"Higher alcohol {ak}",
        "unit": "g/L",
        "initial_value": 0.0,
    })
pd.DataFrame(ic_rows).to_csv(OUT / "initial_conditions.csv", index=False)

# ────────────────────────────────────────────────────────────────────────────
# 9. Kinetic parameters  (kinetic_parameters.csv)
# ────────────────────────────────────────────────────────────────────────────
print("[9/16] Kinetic parameters ...")
_param_set = globals().get("ZENTENO_PARAM_SET", "paper2010")
kinetic_params = {
    "MU0_nom": float(globals().get("MU0_nom", 0.18)),
    "YXN_nom": float(globals().get("YXN_nom", 19.69)),
    "YXG_nom": float(globals().get("YXG_nom", 1.60)),
    "YXF_nom": float(globals().get("YXF_nom", 1.60)),
    "YEG_nom": float(globals().get("YEG_nom", 0.49)),
    "YEF_nom": float(globals().get("YEF_nom", 0.49)),
    "Kn0_nom": float(globals().get("Kn0_nom", 0.01)),
    "Kg0_nom": float(globals().get("Kg0_nom", 7.5)),
    "Kf0_nom": float(globals().get("Kf0_nom", 7.5)),
    "Kig0_nom": float(globals().get("Kig0_nom", 55.0)),
    "Kie0_nom": float(globals().get("Kie0_nom", 40.0)),
    "betaG0_nom": float(globals().get("betaG0_nom", 0.225)),
    "betaF0_nom": float(globals().get("betaF0_nom", 0.225)),
    "Kd0_nom": float(globals().get("Kd0_nom", 0.00044)),
    "MRATE_0": float(globals().get("MRATE_0", 0.0)),
    "T_val": float(globals().get("T_val", 293.15)),
    "R_GAS": float(globals().get("R_GAS", 8.314)),
    "MW_GLU": float(globals().get("MW_GLU", 0.180156)),
    "MW_FRU": float(globals().get("MW_FRU", 0.180156)),
    "MW_ETH": float(globals().get("MW_ETH", 0.046070)),
    "MW_N": float(globals().get("MW_N", 0.014007)),
    "MW_O2": float(globals().get("MW_O2", 0.032)),
    "EPS": float(globals().get("EPS", 1e-9)),
}
kp_df = pd.DataFrame([{"parameter": k, "value": v} for k, v in kinetic_params.items()])
kp_df.to_csv(OUT / "kinetic_parameters.csv", index=False)

# Also save parameter set name for reference
with open(OUT / "kinetic_param_set.txt", "w") as f:
    f.write(f"ZENTENO_PARAM_SET={_param_set}\n")

# ────────────────────────────────────────────────────────────────────────────
# 10. ODE state definitions  (ode_state_definitions.csv)
# ────────────────────────────────────────────────────────────────────────────
print("[10/16] ODE state definitions ...")
ode_rows = [
    {"idx": 1, "name": "X",    "unit": "gDW/L",  "ode_rhs": "v_obj * X  (BIOMASS); 0 (TURNOVER)"},
    {"idx": 2, "name": "N_free","unit": "gN/L",   "ode_rhs": "-vn_eff * X  (BIOMASS); 0 (TURNOVER)"},
    {"idx": 3, "name": "G",    "unit": "g/L",     "ode_rhs": "-MW_GLU * max(0,-v_glu) * X"},
    {"idx": 4, "name": "F",    "unit": "g/L",     "ode_rhs": "-MW_FRU * max(0,-v_fru) * X"},
    {"idx": 5, "name": "E",    "unit": "g/L",     "ode_rhs": "+MW_ETH * max(0, v_eth) * X"},
    {"idx": 6, "name": "O2",   "unit": "g/L",     "ode_rhs": "-MW_O2 * max(0,-v_o2) * X"},
    {"idx": 7, "name": "Prot", "unit": "g/L",     "ode_rhs": "v_obj*X*PROT0 - Prot*Kd + Xa*v_prot - TURNOVER_LAMBDA*Prot"},
    {"idx": 8, "name": "Carb", "unit": "g/L",     "ode_rhs": "v_obj*X*CARB0 - Carb*Kd"},
]
for j, ak in enumerate(_aa_keys):
    ode_rows.append({
        "idx": 9 + j, "name": f"AA_{ak}", "unit": "mmol/L",
        "ode_rhs": "X * v_aa_i  (BIOMASS); 0 (TURNOVER)"
    })
for j, ak in enumerate(_aroma_keys):
    ode_rows.append({
        "idx": 9 + _n_aa + j, "name": f"Aroma_{ak}", "unit": "g/L",
        "ode_rhs": "MW_aroma * v_aroma * X"
    })
pd.DataFrame(ode_rows).to_csv(OUT / "ode_state_definitions.csv", index=False)

# ────────────────────────────────────────────────────────────────────────────
# 11. GAM parameters  (gam_parameters.csv)
# ────────────────────────────────────────────────────────────────────────────
print("[11/16] GAM parameters ...")
gam_rows = {
    "GAM_BASE": float(globals().get("GAM_BASE", 30.49)),
    "GAM_COEFF_P": float(globals().get("GAM_COEFF_P", 16.965)),
    "GAM_COEFF_R": float(globals().get("GAM_COEFF_R", 1.638)),
    "GAM_COEFF_C": float(globals().get("GAM_COEFF_C", 5.210)),
    "Pbase_global": float(globals().get("Pbase_global", 0.46)),
    "Cbase_global": float(globals().get("Cbase_global", 0.37)),
    "Rbase_global": float(globals().get("Rbase_global", 0.06)),
    "PROT_CONTENT_0": float(globals().get("PROT_CONTENT_0", 0.46)),
    "CARB_CONTENT_0": float(globals().get("CARB_CONTENT_0", 0.37)),
    "RNA_FRAC": float(globals().get("RNA_FRAC", 0.06)),
}
pd.DataFrame([{"parameter": k, "value": v} for k, v in gam_rows.items()]).to_csv(
    OUT / "gam_parameters.csv", index=False)

# ────────────────────────────────────────────────────────────────────────────
# 12. Nitrogen config  (nitrogen_config.csv)
# ────────────────────────────────────────────────────────────────────────────
print("[12/16] Nitrogen configuration ...")
_KIN_N = globals().get("KINETIC_N_SOURCE_IDS", [])
_N_atoms = globals().get("N_atoms_map", {})
_N_frac = globals().get("N_frac_map", {})
n_rows = []
for rid in _KIN_N:
    rname = model.reactions.get_by_id(rid).name if rid in rxn_id_to_idx else rid
    n_rows.append({
        "reaction_id": rid,
        "reaction_name": rname,
        "index_1based": rxn_id_to_idx.get(rid, -1),
        "N_atoms": float(_N_atoms.get(rid, 1.0)),
        "N_frac_share": float(_N_frac.get(rid, 0.0)),
    })
pd.DataFrame(n_rows).to_csv(OUT / "nitrogen_config.csv", index=False)

# ────────────────────────────────────────────────────────────────────────────
# 13. Aroma config  (aroma_config.csv)
# ────────────────────────────────────────────────────────────────────────────
print("[13/16] Aroma configuration ...")
_fa = globals().get("found_aromas", {}) or {}
_amw = globals().get("AROMA_MW", {}) or {}
aroma_rows = []
for ak, rid in _fa.items():
    rname = model.reactions.get_by_id(rid).name if rid in rxn_id_to_idx else rid
    aroma_rows.append({
        "aroma_key": ak,
        "reaction_id": rid,
        "reaction_name": rname,
        "index_1based": rxn_id_to_idx.get(rid, -1),
        "MW_g_per_mmol": float(_amw.get(ak, 0.0)),
    })
pd.DataFrame(aroma_rows).to_csv(OUT / "aroma_config.csv", index=False)

# ────────────────────────────────────────────────────────────────────────────
# 14. AA precursor config  (aa_precursor_config.csv)
# ────────────────────────────────────────────────────────────────────────────
print("[14/16] AA precursor configuration ...")
_fp = globals().get("found_precursors", {}) or {}
_aamw = globals().get("AA_MW", {}) or {}
_aalpha = globals().get("AA_ALPHA", {}) or {}
aa_rows = []
for ak, rid in _fp.items():
    rname = model.reactions.get_by_id(rid).name if rid in rxn_id_to_idx else rid
    aa_rows.append({
        "aa_key": ak,
        "reaction_id": rid,
        "reaction_name": rname,
        "index_1based": rxn_id_to_idx.get(rid, -1),
        "MW_g_per_mmol": float(_aamw.get(ak, 0.0)),
        "AA_ALPHA": float(_aalpha.get(ak, 0.0)),
    })
pd.DataFrame(aa_rows).to_csv(OUT / "aa_precursor_config.csv", index=False)

# ────────────────────────────────────────────────────────────────────────────
# 15. Turnover parameters  (turnover_parameters.csv)
# ────────────────────────────────────────────────────────────────────────────
print("[15/16] Turnover parameters ...")
turn_params = {
    "TURNOVER_LAMBDA": float(globals().get("TURNOVER_LAMBDA", 0.03)),
    "K_DEATH": float(globals().get("K_DEATH", 0.005)),
    "VPROT_FLOOR_FRAC": float(globals().get("VPROT_FLOOR_FRAC", 0.90)),
    "XA_FRACTION": float(globals().get("XA_FRACTION", 1.0)),
    "ATPM_LB_NO_GROWTH": float(globals().get("ATPM_LB_NO_GROWTH", 0.7)),
    "ATPM_UB_NO_GROWTH": float(globals().get("ATPM_UB_NO_GROWTH", 1000.0)),
    "BIOMASS_LOCK_FRAC": float(globals().get("BIOMASS_LOCK_FRAC", 0.999)),
    "K_AA_UPTAKE_GROWTH": float(globals().get("K_AA_UPTAKE_GROWTH", 0.08)),
    "AA_UPTAKE_WEIGHT": float(globals().get("AA_UPTAKE_WEIGHT", 1.0)),
    "N_AMMONIA_FRACTION": float(globals().get("N_AMMONIA_FRACTION", 0.50)),
    "N_TOTAL_DEPLETION_THRESHOLD": float(globals().get("N_TOTAL_DEPLETION_THRESHOLD", 1e-3)),
    "APPLY_PRODUCT_CAPS": 1.0 if globals().get("APPLY_PRODUCT_CAPS", False) else 0.0,
}
pd.DataFrame([{"parameter": k, "value": v} for k, v in turn_params.items()]).to_csv(
    OUT / "turnover_parameters.csv", index=False)

# ────────────────────────────────────────────────────────────────────────────
# 16. Tracked bounds log  (tracked_bounds_log.csv)
# ────────────────────────────────────────────────────────────────────────────
print("[16/16] Tracked bounds log ...")
_tracked = globals().get("_tracked_bounds_log", None) or globals().get("tracked_bounds_log", None)
if _tracked and isinstance(_tracked, list):
    pd.DataFrame(_tracked).to_csv(OUT / "tracked_bounds_log.csv", index=False)
    print(f"  tracked_bounds_log.csv: {len(_tracked)} entries")
else:
    # Reconstruct from known set_bounds_logged calls
    log_entries = []
    
    # Anaerobiosis
    log_entries.append({"reaction_id": "r_1992", "lb": 0.0, "ub": 0.0, "group": "Anaerobiosis", "note": "O2 exchange = 0 (anaerobic)"})
    
    # Anaerobic supplements
    for rid in ["r_1757", "r_1915", "r_1994", "r_2106", "r_2134", "r_2137", "r_2189"]:
        rname = model.reactions.get_by_id(rid).name if rid in rxn_id_to_idx else rid
        log_entries.append({"reaction_id": rid, "lb": -1000.0, "ub": 0.0, "group": "Suplementos anaerobicos", "note": rname})
    
    # Carbon sources
    log_entries.append({"reaction_id": "r_1714", "lb": -1000.0, "ub": 0.0, "group": "Fuentes de carbono", "note": "Glucosa solo uptake"})
    log_entries.append({"reaction_id": "r_1709", "lb": -1000.0, "ub": 0.0, "group": "Fuentes de carbono", "note": "Fructosa solo uptake"})
    
    # ATPM
    log_entries.append({"reaction_id": "r_4046", "lb": 0.7, "ub": 0.7, "group": "ATPM", "note": "ATPM fijo (base model)"})
    
    # Products
    for rid in ["r_1761", "r_1808", "r_1634", "r_2056"]:
        rname = model.reactions.get_by_id(rid).name if rid in rxn_id_to_idx else rid
        log_entries.append({"reaction_id": rid, "lb": 0.0, "ub": 1000.0, "group": "Productos abiertos", "note": rname})
    
    # Water/protons
    log_entries.append({"reaction_id": "r_2100", "lb": -1000.0, "ub": 1000.0, "group": "Agua y protones", "note": "water exchange"})
    log_entries.append({"reaction_id": "r_1832", "lb": -1000.0, "ub": 1000.0, "group": "Agua y protones", "note": "H+ exchange"})
    
    # Leak closure
    for rid, note in [("r_4501","turanose"), ("r_1650","trehalose"), ("r_4502","G6P"), ("r_4538","6PG")]:
        log_entries.append({"reaction_id": rid, "lb": 0.0, "ub": 0.0, "group": "Leak closure", "note": note})
    
    # Alt sink caps
    for rid, cap in [("r_1549",0.5), ("r_4659",0.1), ("r_4731",0.1), ("r_4737",0.1), ("r_1989",0.0)]:
        rname = model.reactions.get_by_id(rid).name if rid in rxn_id_to_idx else rid
        log_entries.append({"reaction_id": rid, "lb": 0.0, "ub": cap, "group": "Alt sink caps", "note": rname})
    
    pd.DataFrame(log_entries).to_csv(OUT / "tracked_bounds_log.csv", index=False)
    print(f"  tracked_bounds_log.csv (reconstructed): {len(log_entries)} entries")

# ────────────────────────────────────────────────────────────────────────────
# Summary
# ────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 80)
print("[dFVB Export] COMPLETE")
print("=" * 80)
exported_files = sorted(OUT.glob("*.csv"))
for f in exported_files:
    sz = f.stat().st_size
    print(f"  {f.name:45s}  ({sz:>8,d} bytes)")
print(f"\nTotal: {len(exported_files)} files in {OUT.resolve()}")
print("Nota: Para importar S_matrix en MATLAB, usar S_matrix_sparse.csv (formato triplet)")
