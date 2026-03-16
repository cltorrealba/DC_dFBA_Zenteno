# Notebook 1 — snippet export simultáneo
# Ejecutar DESPUÉS de la celda final de GAM variable / pairwise / ethyl acetate
# Produce: out/dfba_vargam_metadata.jl
#
# Extiende el snippet anterior agregando:
#   - baseline_time_h
#   - baseline_gam_mmol_gdw
# para que Notebook 2 pueda construir GAM_EXTRA_FE.

from pathlib import Path
import numpy as np

required = [
    "OUT_DIR",
    "OBJ_ID", "GLU_ID", "FRU_ID", "ETH_ID", "O2_ID", "ATPM_ID",
    "PROT_RXN_ID",
    "KINETIC_N_SOURCE_IDS", "N_VITAMIN_EX_IDS",
    "found_precursors", "found_aromas",
    "N_atoms_map", "N_frac_map",
    "AA_ALPHA", "AA_MW", "AROMA_MW",
    "PROT_CONTENT_0", "CARB_CONTENT_0", "RNA_FRAC",
    "K_DEATH", "TURNOVER_LAMBDA", "XA_FRACTION",
    "N_AMMONIA_FRACTION", "N_TOTAL_DEPLETION_THRESHOLD",
    "K_AA_UPTAKE_GROWTH",
    "ATPM_LB_NO_GROWTH", "ATPM_UB_NO_GROWTH",
    "GAM_BASE", "GAM_COEFF_P", "GAM_COEFF_R", "GAM_COEFF_C",
    "Pbase_global", "Cbase_global", "Rbase_global",
    "_find_biomass_rxn",
]
missing = [k for k in required if k not in globals()]
if missing:
    raise RuntimeError(f"Faltan símbolos requeridos para exportar metadata: {missing}")

META_FILE = Path(OUT_DIR) / "dfba_vargam_metadata.jl"

def to_julia(obj):
    if isinstance(obj, str):
        return '"' + obj.replace("\\", "\\\\").replace('"', '\\"') + '"'
    if obj is None:
        return "nothing"
    if isinstance(obj, bool):
        return "true" if obj else "false"
    if isinstance(obj, (int, np.integer)):
        return str(int(obj))
    if isinstance(obj, (float, np.floating)):
        v = float(obj)
        if np.isnan(v):
            return "NaN"
        if np.isinf(v):
            return "Inf" if v > 0 else "-Inf"
        return repr(v)
    if isinstance(obj, dict):
        return "Dict(" + ", ".join(f'{to_julia(str(k))} => {to_julia(v)}' for k, v in obj.items()) + ")"
    if isinstance(obj, (list, tuple)):
        return "[" + ", ".join(to_julia(v) for v in obj) + "]"
    return to_julia(str(obj))

pairwise_constraints = []
if "active_pairs" in globals() and active_pairs is not None and len(active_pairs) > 0:
    for _, row in active_pairs.iterrows():
        pairwise_constraints.append((
            str(row["ester_rid"]),
            str(row["alcohol_rid"]),
            float(row["phi_pair"]),
        ))

ethyl_acetate_soft = None
if globals().get("INCLUDE_ETHYL_ACETATE_SOFT_COUPLING", False):
    ethyl_acetate_soft = {
        "ester_rid": str(globals().get("ETHYL_ACETATE_RID_EA", "r_1765")),
        "alcohol_rid": str(globals().get("ETHANOL_RID_EA", ETH_ID)),
        "phi": float(globals().get("PHI_ETHYL_ACETATE_STATIC", 0.005)),
    }

bio_rxn = _find_biomass_rxn(model)
targets = {"atp", "adp", "h2o", "h+", "phosphate", "water", "proton"}
gam_targets = []
for met, coef in bio_rxn.metabolites.items():
    nm = met.name.lower()
    if any(t in nm for t in targets):
        gam_targets.append((met.id, float(np.sign(coef))))

baseline_time_h = []
baseline_gam_mmol_gdw = []

if "states_df_phase" in globals() and states_df_phase is not None and len(states_df_phase) > 0:
    if "t_h" in states_df_phase.columns and "GAM_mmol_gDW" in states_df_phase.columns:
        baseline_time_h = [float(x) for x in states_df_phase["t_h"].tolist()]
        baseline_gam_mmol_gdw = [float(x) for x in states_df_phase["GAM_mmol_gDW"].tolist()]
elif "states_df" in globals() and states_df is not None and len(states_df) > 0:
    if "t_h" in states_df.columns and "GAM_mmol_gDW" in states_df.columns:
        baseline_time_h = [float(x) for x in states_df["t_h"].tolist()]
        baseline_gam_mmol_gdw = [float(x) for x in states_df["GAM_mmol_gDW"].tolist()]

meta = {
    "obj_id": str(OBJ_ID),
    "glu_id": str(GLU_ID),
    "fru_id": str(FRU_ID),
    "eth_id": str(ETH_ID),
    "o2_id": str(O2_ID),
    "atpm_id": str(ATPM_ID),
    "prot_rxn_id": str(PROT_RXN_ID),
    "kinetic_n_source_ids": [str(x) for x in KINETIC_N_SOURCE_IDS],
    "vitamin_n_source_ids": [str(x) for x in N_VITAMIN_EX_IDS],
    "aa_exchange_ids": {str(k): str(v) for k, v in found_precursors.items()},
    "aroma_exchange_ids": {str(k): str(v) for k, v in found_aromas.items()},
    "n_atoms_map": {str(k): float(v) for k, v in N_atoms_map.items()},
    "n_frac_map": {str(k): float(v) for k, v in N_frac_map.items()},
    "aa_alpha": {str(k): float(v) for k, v in AA_ALPHA.items()},
    "aa_mw": {str(k): float(v) for k, v in AA_MW.items()},
    "aroma_mw": {str(k): float(v) for k, v in AROMA_MW.items()},
    "PROT_CONTENT_0": float(PROT_CONTENT_0),
    "CARB_CONTENT_0": float(CARB_CONTENT_0),
    "RNA_FRAC": float(RNA_FRAC),
    "K_DEATH": float(K_DEATH),
    "TURNOVER_LAMBDA": float(TURNOVER_LAMBDA),
    "XA_FRACTION": float(XA_FRACTION),
    "N_AMMONIA_FRACTION": float(N_AMMONIA_FRACTION),
    "N_TOTAL_DEPLETION_THRESHOLD": float(N_TOTAL_DEPLETION_THRESHOLD),
    "K_AA_UPTAKE_GROWTH": float(K_AA_UPTAKE_GROWTH),
    "ATPM_LB_NO_GROWTH": float(ATPM_LB_NO_GROWTH),
    "ATPM_UB_NO_GROWTH": float(ATPM_UB_NO_GROWTH),
    "GAM_BASE": float(GAM_BASE),
    "GAM_COEFF_P": float(GAM_COEFF_P),
    "GAM_COEFF_R": float(GAM_COEFF_R),
    "GAM_COEFF_C": float(GAM_COEFF_C),
    "Pbase_global": float(Pbase_global),
    "Cbase_global": float(Cbase_global),
    "Rbase_global": float(Rbase_global),
    "pairwise_constraints": pairwise_constraints,
    "baseline_time_h": baseline_time_h,
    "baseline_gam_mmol_gdw": baseline_gam_mmol_gdw,
    "gam_callback": {
        "biomass_rxn_id": str(bio_rxn.id),
        "targets": gam_targets,
    },
}
if ethyl_acetate_soft is not None:
    meta["ethyl_acetate_soft"] = ethyl_acetate_soft

content = "const DFBA_META = " + to_julia(meta) + "\n"
META_FILE.write_text(content, encoding="utf-8")

print(f"Metadata simultánea exportada en: {META_FILE}")
print("Pairwise:", meta["pairwise_constraints"])
print("Baseline GAM points:", len(meta["baseline_time_h"]))
