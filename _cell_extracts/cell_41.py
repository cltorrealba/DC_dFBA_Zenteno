import pandas as pd
import numpy as np
import cobra
import math
import re

# ======================================================================================
# Acetate esters: REFORMULACIÓN PAIRWISE (sin dominancia de ethyl acetate)
# Análisis snapshot-level sobre la trayectoria baseline – NO es re-simulación dinámica
# ======================================================================================

print("\n" + "=" * 100)
print("Acetate esters reformulados: restricciones pairwise (sin ethanol global)")
print("  [NOTA: snapshot-level sobre trayectoria baseline, no re-simulación dinámica completa]")
print("=" * 100)

required_symbols_pairwise = [
	"model", "OUT_DIR",
	"states_df_phase", "flux_df_phase", "phase_meta_df",
	"_build_phase_snapshot_model",
	"available_aroma_targets",
	"OBJ_ID", "ATPM_ID", "ETH_ID",
	"acetate_ester_targets",
]
missing_symbols_pairwise = [s for s in required_symbols_pairwise if s not in globals()]
if missing_symbols_pairwise:
	print(f"[WARN] Símbolos faltantes: {missing_symbols_pairwise}")

# ─────────────────────────────────────────────────────────────────────────────────────
# CONFIGURACIÓN
# ─────────────────────────────────────────────────────────────────────────────────────

ENABLE_PAIRWISE_COUPLING      = True
USE_PAIRWISE_ESTER_COUPLING   = True
EXCLUDE_ETHANOL_FROM_COUPLING = True        # No incluir ETH_ID en acoplamiento de esters secundarios
INCLUDE_ETHYL_ACETATE_PAIR    = False       # ethyl acetate sin restricción pairwise dedicada
ETHYL_ACETATE_RID             = "r_1765"   # ethyl acetate exchange (mantener fuera del acoplamiento)
PAIRWISE_PHI_MAX              = 0.15
PAIRWISE_PHI_GRID             = [0.00, 0.03, 0.06, 0.12, 0.15]
PAIRWISE_PHI_STATIC           = 0.08       # φ fijo para cada par

ZERO_TOL       = 1e-10   # tolerancia numérica para limpiar ruido
ESTER_BIND_TOL = 1e-7    # tolerancia para declarar binding

RUN_PAIRWISE_SENSITIVITY = True
SAVE_PAIRWISE_CSV        = True

print(f"USE_PAIRWISE_ESTER_COUPLING    = {USE_PAIRWISE_ESTER_COUPLING}")
print(f"EXCLUDE_ETHANOL_FROM_COUPLING  = {EXCLUDE_ETHANOL_FROM_COUPLING}")
print(f"INCLUDE_ETHYL_ACETATE_PAIR     = {INCLUDE_ETHYL_ACETATE_PAIR}")
print(f"PAIRWISE_PHI_MAX               = {PAIRWISE_PHI_MAX}")
print(f"PAIRWISE_PHI_STATIC            = {PAIRWISE_PHI_STATIC}")
print(f"ZERO_TOL                       = {ZERO_TOL}")
print(f"ESTER_BIND_TOL                 = {ESTER_BIND_TOL}")

# ─────────────────────────────────────────────────────────────────────────────────────
# MAPEO MANUAL CORREGIDO: ester_rid -> alcohol_rid
# IDs verificados manualmente en yeast-GEM
#   r_1862 = isoamyl acetate exchange    -> r_1865 = isoamylol exchange
#   r_1867 = isobutyl acetate exchange   -> r_1866 = isobutanol exchange
#   r_2000 = phenethyl acetate exchange  -> r_1589 = 2-phenylethanol exchange
# ethyl acetate (r_1765) se mantiene fuera del acoplamiento pairwise principal
# ─────────────────────────────────────────────────────────────────────────────────────

manual_alcohol_map = {
	"r_1862": "r_1865",   # isoamyl acetate  -> isoamylol
	"r_1867": "r_1866",   # isobutyl acetate -> isobutanol
	"r_2000": "r_1589",   # phenethyl acetate -> 2-phenylethanol
}

# Etiquetas para FVA y verificación de signo
PAIR_LABELS = {
	"r_1862": "isoamyl acetate",
	"r_1865": "isoamylol",
	"r_1867": "isobutyl acetate",
	"r_1866": "isobutanol",
	"r_2000": "phenethyl acetate",
	"r_1589": "2-phenylethanol",
}
ALCOHOL_SIGN_IDS = ["r_1865", "r_1866", "r_1589"]

print(f"\nMapeo manual corregido (ester_rid -> alcohol_rid):")
for k_rid, v_rid in manual_alcohol_map.items():
	print(f"  {k_rid} ({PAIR_LABELS[k_rid]}) -> {v_rid} ({PAIR_LABELS[v_rid]})")

# ─────────────────────────────────────────────────────────────────────────────────────
# A. RECONSTRUIR MAPEO ESTER ↔ ALCOHOL (excluyendo ethanol)
# ─────────────────────────────────────────────────────────────────────────────────────

rxn_ids_all = {r.id for r in model.reactions}
ester_alcohol_pairwise_rows = []

for ester_rid, ester_label in acetate_ester_targets.items():
	# Resolución: primero manual, luego fallback a ester_to_alcohol si existe
	alcohol_rid = manual_alcohol_map.get(ester_rid)
	if alcohol_rid is None and "ester_to_alcohol" in globals():
		alcohol_rid = ester_to_alcohol.get(ester_rid)

	if alcohol_rid is None:
		ester_alcohol_pairwise_rows.append({
			"ester_rid": ester_rid, "ester_label": ester_label,
			"alcohol_rid": None, "alcohol_name": None,
			"pair_active": False, "exclusion_reason": "no_alcohol_found", "phi_pair": 0.0,
		})
		continue

	# Filtro 1: excluir ethanol del acoplamiento principal
	is_ethanol_pair = (alcohol_rid == ETH_ID) or ("ethanol" in str(ester_label).lower())
	if EXCLUDE_ETHANOL_FROM_COUPLING and is_ethanol_pair:
		alcohol_name = model.reactions.get_by_id(alcohol_rid).name if alcohol_rid in rxn_ids_all else None
		ester_alcohol_pairwise_rows.append({
			"ester_rid": ester_rid, "ester_label": ester_label,
			"alcohol_rid": alcohol_rid, "alcohol_name": alcohol_name,
			"pair_active": False, "exclusion_reason": "ethanol_excluded", "phi_pair": 0.0,
		})
		continue

	# Filtro 2: excluir ethyl acetate
	if not INCLUDE_ETHYL_ACETATE_PAIR and ester_rid == ETHYL_ACETATE_RID:
		alcohol_name = model.reactions.get_by_id(alcohol_rid).name if alcohol_rid in rxn_ids_all else None
		ester_alcohol_pairwise_rows.append({
			"ester_rid": ester_rid, "ester_label": ester_label,
			"alcohol_rid": alcohol_rid, "alcohol_name": alcohol_name,
			"pair_active": False, "exclusion_reason": "ethyl_acetate_excluded", "phi_pair": 0.0,
		})
		continue

	# Filtro 3: verificar existencia en modelo
	if ester_rid not in rxn_ids_all or alcohol_rid not in rxn_ids_all:
		ester_alcohol_pairwise_rows.append({
			"ester_rid": ester_rid, "ester_label": ester_label,
			"alcohol_rid": alcohol_rid, "alcohol_name": None,
			"pair_active": False, "exclusion_reason": "missing_from_model", "phi_pair": 0.0,
		})
		continue

	# Par activo
	ester_alcohol_pairwise_rows.append({
		"ester_rid": ester_rid, "ester_label": ester_label,
		"alcohol_rid": alcohol_rid,
		"alcohol_name": model.reactions.get_by_id(alcohol_rid).name,
		"pair_active": True, "exclusion_reason": None,
		"phi_pair": float(PAIRWISE_PHI_STATIC),
	})

ester_alcohol_pairwise_df = pd.DataFrame(ester_alcohol_pairwise_rows)
active_pairs = ester_alcohol_pairwise_df[ester_alcohol_pairwise_df["pair_active"] == True].copy()

print(f"\nTotales:")
print(f"  Pares identificados:   {len(ester_alcohol_pairwise_df)}")
print(f"  Pares activos:         {len(active_pairs)}")
print(f"  Excluidos (ethanol):   {(ester_alcohol_pairwise_df['exclusion_reason'] == 'ethanol_excluded').sum()}")
print(f"  Excluidos (ethylac.):  {(ester_alcohol_pairwise_df['exclusion_reason'] == 'ethyl_acetate_excluded').sum()}")
print(f"  Excluidos (model):     {(ester_alcohol_pairwise_df['exclusion_reason'] == 'missing_from_model').sum()}")
print(f"  Sin alcohol mapeado:   {(ester_alcohol_pairwise_df['exclusion_reason'] == 'no_alcohol_found').sum()}")

print("\n[Mapeo pairwise final]")
display(ester_alcohol_pairwise_df)

active_ester_ids_pairwise   = list(active_pairs["ester_rid"].unique())
active_alcohol_ids_pairwise = list(active_pairs["alcohol_rid"].unique())

# ─────────────────────────────────────────────────────────────────────────────────────
# B. FUNCIONES AUXILIARES PAIRWISE
# ─────────────────────────────────────────────────────────────────────────────────────

def _add_pairwise_ester_coupling_constraints(mtmp, pairwise_df):
	"""
	Agrega restricciones: v_ester >= phi * v_alcohol  para cada par activo.
	  ⟺  v_ester - phi * v_alcohol >= 0

	Retorna (constraints_added, active_map, phi_map)
	  active_map: {(ester_rid, alcohol_rid): constraint_object}
	  phi_map:    {(ester_rid, alcohol_rid): phi_pair}
	"""
	if not USE_PAIRWISE_ESTER_COUPLING or pairwise_df.empty:
		return [], {}, {}

	rb = {r.id: r for r in mtmp.reactions}
	constraints_added = []
	active_map = {}
	phi_map    = {}  # (ester_rid, alcohol_rid) -> phi_pair

	for _, row_p in pairwise_df[pairwise_df["pair_active"] == True].iterrows():
		ester_rid_p  = str(row_p["ester_rid"])
		alcohol_rid_p = str(row_p["alcohol_rid"])
		phi_val_p     = float(row_p["phi_pair"])

		if phi_val_p <= 0.0:
			continue
		if ester_rid_p not in rb or alcohol_rid_p not in rb:
			continue

		ester_expr   = rb[ester_rid_p].flux_expression
		alcohol_expr = rb[alcohol_rid_p].flux_expression

		cname = f"pair_{ester_rid_p}__{alcohol_rid_p}_phi{str(phi_val_p).replace('.', '_')}"
		cons  = mtmp.problem.Constraint(
			ester_expr - phi_val_p * alcohol_expr,
			lb=0.0,
			name=cname
		)
		mtmp.add_cons_vars([cons])
		constraints_added.append(cons)
		active_map[(ester_rid_p, alcohol_rid_p)] = cons
		phi_map[(ester_rid_p, alcohol_rid_p)]    = phi_val_p

	return constraints_added, active_map, phi_map


def _build_phase_snapshot_model_pairwise(t_phase, *, enable_pairwise=True, pairwise_df=None):
	"""
	Wrapper snapshot-level sobre _build_phase_snapshot_model().
	Añade restricciones pairwise al snapshot; NO re-simula la dinámica completa.

	El binding se calcula con el phi_pair correcto de cada par (via phi_map),
	usando: lhs_margin = v_ester - phi_pair * v_alcohol  (binding si <= ESTER_BIND_TOL)
	"""
	snap         = _build_phase_snapshot_model(float(t_phase))
	m_base       = snap["model"]
	sol_base     = snap["solution"]
	state_row    = snap["state_row"]
	mode         = snap["mode"]
	objective_id = snap["objective_id"]

	m_pairwise   = m_base.copy()
	sol_pairwise = sol_base

	constraints_pairwise = []
	active_pair_map      = {}
	phi_map              = {}
	constraints_binding  = {}

	if enable_pairwise and pairwise_df is not None and not pairwise_df.empty:
		constraints_pairwise, active_pair_map, phi_map = \
			_add_pairwise_ester_coupling_constraints(m_pairwise, pairwise_df)

		try:
			sol_pairwise = m_pairwise.optimize()
		except Exception as exc:
			print(f"[WARN] t={t_phase}: error en optimize con pairwise: {exc}")
			sol_pairwise = sol_base

		if getattr(sol_pairwise, "status", "") == "optimal":
			for (ester_rid_b, alcohol_rid_b) in active_pair_map:
				phi_pair_b = phi_map.get((ester_rid_b, alcohol_rid_b), 0.0)

				v_ester_b = float(sol_pairwise.fluxes.get(ester_rid_b, 0.0))
				v_alc_b   = float(sol_pairwise.fluxes.get(alcohol_rid_b, 0.0))

				# limpiar ruido numérico
				if abs(v_ester_b) < ZERO_TOL: v_ester_b = 0.0
				if abs(v_alc_b)   < ZERO_TOL: v_alc_b   = 0.0

				lhs_margin = v_ester_b - phi_pair_b * v_alc_b
				if abs(lhs_margin) < ZERO_TOL: lhs_margin = 0.0

				constraints_binding[(ester_rid_b, alcohol_rid_b)] = {
					"phi_pair":   phi_pair_b,
					"v_ester":    v_ester_b,
					"v_alcohol":  v_alc_b,
					"lhs_margin": lhs_margin,
					"binding":    lhs_margin <= ESTER_BIND_TOL,
				}

	return {
		"t_h":               float(t_phase),
		"mode":              mode,
		"objective_id":      objective_id,
		"state_row":         state_row,
		"model_baseline":    m_base,
		"model_pairwise":    m_pairwise,
		"solution_baseline": sol_base,
		"solution_pairwise": sol_pairwise,
		"constraints_pairwise": constraints_pairwise,
		"active_pair_map":   active_pair_map,
		"phi_map":           phi_map,
		"constraints_binding": constraints_binding,
	}


# ─────────────────────────────────────────────────────────────────────────────────────
# C. COMPARACIÓN BASELINE VS PAIRWISE EN EL TIEMPO
# (snapshot-level sobre trayectoria baseline)
# ─────────────────────────────────────────────────────────────────────────────────────

analysis_times_pw = flux_df_phase["t_h"].to_numpy(dtype=float)
state_times_pw    = states_df_phase["t_h"].to_numpy(dtype=float)

comparison_rows_pw  = []
snapshot_cache_pairwise = {}

for k, t_k in enumerate(analysis_times_pw):
	res_pw = _build_phase_snapshot_model_pairwise(
		t_k,
		enable_pairwise=ENABLE_PAIRWISE_COUPLING,
		pairwise_df=ester_alcohol_pairwise_df,
	)

	snapshot_cache_pairwise[float(t_k)] = {
		"model":    res_pw["model_pairwise"],
		"solution": res_pw["solution_pairwise"],
	}

	row_state = res_pw["state_row"]
	sol_b     = res_pw["solution_baseline"]
	sol_pw    = res_pw["solution_pairwise"]
	obj_id    = res_pw["objective_id"]

	mu_b  = float(sol_b.fluxes.get(OBJ_ID, np.nan))  if getattr(sol_b,  "status", "") == "optimal" else np.nan
	mu_pw = float(sol_pw.fluxes.get(OBJ_ID, np.nan)) if getattr(sol_pw, "status", "") == "optimal" else np.nan
	obj_b  = float(sol_b.fluxes.get(obj_id, np.nan))  if getattr(sol_b,  "status", "") == "optimal" else np.nan
	obj_pw = float(sol_pw.fluxes.get(obj_id, np.nan)) if getattr(sol_pw, "status", "") == "optimal" else np.nan

	out_row = {
		"t_h":                      float(t_k),
		"mode":                     res_pw["mode"],
		"objective_id":             obj_id,
		"objective_value_baseline": obj_b,
		"objective_value_pairwise": obj_pw,
		"delta_objective":          obj_pw - obj_b if np.isfinite(obj_b) and np.isfinite(obj_pw) else np.nan,
		"mu_h_inv_baseline":        mu_b,
		"mu_h_inv_pairwise":        mu_pw,
		"delta_mu":                 mu_pw - mu_b if np.isfinite(mu_b) and np.isfinite(mu_pw) else np.nan,
		"n_constraints_active":     len(res_pw["constraints_pairwise"]),
		"n_constraints_binding":    sum(1 for v in res_pw["constraints_binding"].values() if v["binding"]),
		"N_gN_L":  float(row_state["N_gN_L"])  if "N_gN_L"  in row_state.index else np.nan,
		"E_g_L":   float(row_state["E_g_L"])   if "E_g_L"   in row_state.index else np.nan,
		"O2_g_L":  float(row_state["O2_g_L"])  if "O2_g_L"  in row_state.index else np.nan,
		"X_gDW_L": float(row_state["X_gDW_L"]) if "X_gDW_L" in row_state.index else np.nan,
	}

	# Flujos por éster
	for ester_rid_c, ester_label_c in acetate_ester_targets.items():
		vb  = float(sol_b.fluxes.get(ester_rid_c, np.nan))  if getattr(sol_b,  "status", "") == "optimal" else np.nan
		vpw = float(sol_pw.fluxes.get(ester_rid_c, np.nan)) if getattr(sol_pw, "status", "") == "optimal" else np.nan
		if np.isfinite(vb)  and abs(vb)  < ZERO_TOL: vb  = 0.0
		if np.isfinite(vpw) and abs(vpw) < ZERO_TOL: vpw = 0.0
		out_row[f"ester_{ester_rid_c}_baseline"] = vb
		out_row[f"ester_{ester_rid_c}_pairwise"] = vpw
		out_row[f"ester_{ester_rid_c}_delta"]    = vpw - vb if np.isfinite(vb) and np.isfinite(vpw) else np.nan

	# Flujos por alcohol
	for alcohol_rid_c in active_alcohol_ids_pairwise:
		vb_alc  = float(sol_b.fluxes.get(alcohol_rid_c, np.nan))  if getattr(sol_b,  "status", "") == "optimal" else np.nan
		vpw_alc = float(sol_pw.fluxes.get(alcohol_rid_c, np.nan)) if getattr(sol_pw, "status", "") == "optimal" else np.nan
		if np.isfinite(vb_alc)  and abs(vb_alc)  < ZERO_TOL: vb_alc  = 0.0
		if np.isfinite(vpw_alc) and abs(vpw_alc) < ZERO_TOL: vpw_alc = 0.0
		out_row[f"alcohol_{alcohol_rid_c}_baseline"] = vb_alc
		out_row[f"alcohol_{alcohol_rid_c}_pairwise"] = vpw_alc
		out_row[f"alcohol_{alcohol_rid_c}_delta"]    = vpw_alc - vb_alc if np.isfinite(vb_alc) and np.isfinite(vpw_alc) else np.nan

	# Binding y lhs_margin por par
	for (est_id, alc_id), binfo in res_pw["constraints_binding"].items():
		out_row[f"binding_{est_id}__{alc_id}"] = binfo["binding"]
		out_row[f"lhs_{est_id}__{alc_id}"]     = binfo["lhs_margin"]

	comparison_rows_pw.append(out_row)

ester_pairwise_time_df = pd.DataFrame(comparison_rows_pw)

print("\n[Comparativa por tiempo: baseline vs pairwise]")
display(ester_pairwise_time_df[[
	"t_h", "mode", "objective_value_baseline", "objective_value_pairwise",
	"delta_objective", "mu_h_inv_baseline", "mu_h_inv_pairwise",
	"n_constraints_active", "n_constraints_binding"
]].head(10))

# ─────────────────────────────────────────────────────────────────────────────────────
# D. DIAGNÓSTICO EXPLÍCITO POR PAR
# Snapshot-level sobre trayectoria baseline; columnas completas requeridas
# ─────────────────────────────────────────────────────────────────────────────────────

ester_alcohol_rows = []

for k_d, t_k_d in enumerate(analysis_times_pw):
	tk_rows = ester_pairwise_time_df[ester_pairwise_time_df["t_h"] == float(t_k_d)]
	if tk_rows.empty:
		continue
	tr_d   = tk_rows.iloc[0]
	mode_d = tr_d["mode"]

	for _, pr in active_pairs.iterrows():
		ester_rid_d   = str(pr["ester_rid"])
		ester_label_d = str(pr["ester_label"])
		alcohol_rid_d = str(pr["alcohol_rid"])
		alcohol_nm_d  = str(pr["alcohol_name"])
		phi_pair_d    = float(pr["phi_pair"])

		vb_est  = float(tr_d.get(f"ester_{ester_rid_d}_baseline",    np.nan))
		vpw_est = float(tr_d.get(f"ester_{ester_rid_d}_pairwise",    np.nan))
		vb_alc  = float(tr_d.get(f"alcohol_{alcohol_rid_d}_baseline", np.nan))
		vpw_alc = float(tr_d.get(f"alcohol_{alcohol_rid_d}_pairwise", np.nan))

		# limpiar ruido numérico
		if np.isfinite(vb_est)  and abs(vb_est)  < ZERO_TOL: vb_est  = 0.0
		if np.isfinite(vpw_est) and abs(vpw_est) < ZERO_TOL: vpw_est = 0.0
		if np.isfinite(vb_alc)  and abs(vb_alc)  < ZERO_TOL: vb_alc  = 0.0
		if np.isfinite(vpw_alc) and abs(vpw_alc) < ZERO_TOL: vpw_alc = 0.0

		ratio_baseline = (vb_est  / vb_alc)  if (np.isfinite(vb_alc)  and abs(vb_alc)  > ZERO_TOL) else np.nan
		ratio_pairwise = (vpw_est / vpw_alc) if (np.isfinite(vpw_alc) and abs(vpw_alc) > ZERO_TOL) else np.nan

		lhs_margin = vpw_est - phi_pair_d * vpw_alc if (np.isfinite(vpw_est) and np.isfinite(vpw_alc)) else np.nan
		if np.isfinite(lhs_margin) and abs(lhs_margin) < ZERO_TOL:
			lhs_margin = 0.0
		binding_d = (lhs_margin <= ESTER_BIND_TOL) if np.isfinite(lhs_margin) else False

		ester_alcohol_rows.append({
			"t_h":                float(t_k_d),
			"mode":               mode_d,
			"ester_rid":          ester_rid_d,
			"ester_label":        ester_label_d,
			"alcohol_rid":        alcohol_rid_d,
			"alcohol_name":       alcohol_nm_d,
			"phi_pair":           phi_pair_d,
			"v_alcohol_baseline": vb_alc,
			"v_alcohol_pairwise": vpw_alc,
			"v_ester_baseline":   vb_est,
			"v_ester_pairwise":   vpw_est,
			"ratio_baseline":     ratio_baseline,
			"ratio_pairwise":     ratio_pairwise,
			"lhs_margin":         lhs_margin,
			"binding":            binding_d,
		})

ester_pairwise_pair_diagnostics_df = pd.DataFrame(ester_alcohol_rows)

print("\n[Diagnóstico explícito por par (snapshot-level sobre trayectoria baseline)]")
diag_cols = [
	"t_h", "mode", "ester_rid", "ester_label",
	"alcohol_rid", "alcohol_name", "phi_pair",
	"v_alcohol_baseline", "v_alcohol_pairwise",
	"v_ester_baseline",   "v_ester_pairwise",
	"ratio_baseline", "ratio_pairwise",
	"lhs_margin", "binding",
]
display(ester_pairwise_pair_diagnostics_df[diag_cols])

# ─────────────────────────────────────────────────────────────────────────────────────
# E. FVA BASELINE VS PAIRWISE – 3 pares activos + verificación de signo en alcoholes
# Ejecutado solo en tiempos de fase (phase_meta_df); fraction_of_optimum = 0.95
# ─────────────────────────────────────────────────────────────────────────────────────

from cobra.flux_analysis import flux_variability_analysis as _fva

FVA_ALL_PAIR_IDS = list(PAIR_LABELS.keys())   # 6 reactions: 3 esters + 3 alcohols

phase_times_fva = sorted(
	set(float(x) for x in phase_meta_df["t_h"].tolist()) if not phase_meta_df.empty else set()
)
print(f"\nTiempos de fase para FVA: {phase_times_fva}")

fva_pair_rows = []

for t_fva in phase_times_fva:
	# Obtener modelos de los caches ya construidos
	snap_b_fva  = snapshot_cache_baseline.get(float(t_fva), {}) if "snapshot_cache_baseline" in globals() else {}
	snap_pw_fva = snapshot_cache_pairwise.get(float(t_fva), {})

	m_b_fva  = snap_b_fva.get("model")
	m_pw_fva = snap_pw_fva.get("model")

	if m_b_fva is None or m_pw_fva is None:
		print(f"[WARN] t={t_fva}: modelo no encontrado en cache (baseline={m_b_fva is not None}, pairwise={m_pw_fva is not None})")
		continue

	# Solo las IDs que existen en ambos modelos
	rb_b  = {r.id for r in m_b_fva.reactions}
	rb_pw = {r.id for r in m_pw_fva.reactions}
	fva_ids_fva = [rid for rid in FVA_ALL_PAIR_IDS if rid in rb_b and rid in rb_pw]

	if not fva_ids_fva:
		print(f"[WARN] t={t_fva}: ninguna reacción del par encontrada en modelo.")
		continue

	mode_fva = ""
	if not phase_meta_df.empty and "mode" in phase_meta_df.columns:
		_rows = phase_meta_df.loc[phase_meta_df["t_h"] == float(t_fva), "mode"]
		if not _rows.empty:
			mode_fva = str(_rows.iloc[0])

	try:
		fva_b_res  = _fva(m_b_fva,  reaction_list=fva_ids_fva, fraction_of_optimum=0.95)
	except Exception as exc:
		print(f"[WARN] FVA baseline t={t_fva}: {exc}")
		fva_b_res  = pd.DataFrame()

	try:
		fva_pw_res = _fva(m_pw_fva, reaction_list=fva_ids_fva, fraction_of_optimum=0.95)
	except Exception as exc:
		print(f"[WARN] FVA pairwise t={t_fva}: {exc}")
		fva_pw_res = pd.DataFrame()

	for rid in fva_ids_fva:
		fva_min_b  = float(fva_b_res.loc[rid, "minimum"])  if not fva_b_res.empty  and rid in fva_b_res.index  else np.nan
		fva_max_b  = float(fva_b_res.loc[rid, "maximum"])  if not fva_b_res.empty  and rid in fva_b_res.index  else np.nan
		fva_min_pw = float(fva_pw_res.loc[rid, "minimum"]) if not fva_pw_res.empty and rid in fva_pw_res.index else np.nan
		fva_max_pw = float(fva_pw_res.loc[rid, "maximum"]) if not fva_pw_res.empty and rid in fva_pw_res.index else np.nan

		# limpiar ruido numérico
		for _v, _nm in [(fva_min_b, "fva_min_b"), (fva_max_b, "fva_max_b"),
		                (fva_min_pw, "fva_min_pw"), (fva_max_pw, "fva_max_pw")]:
			pass  # asignamos abajo
		if np.isfinite(fva_min_b)  and abs(fva_min_b)  < ZERO_TOL: fva_min_b  = 0.0
		if np.isfinite(fva_max_b)  and abs(fva_max_b)  < ZERO_TOL: fva_max_b  = 0.0
		if np.isfinite(fva_min_pw) and abs(fva_min_pw) < ZERO_TOL: fva_min_pw = 0.0
		if np.isfinite(fva_max_pw) and abs(fva_max_pw) < ZERO_TOL: fva_max_pw = 0.0

		range_b  = fva_max_b  - fva_min_b  if all(np.isfinite([fva_min_b,  fva_max_b]))  else np.nan
		range_pw = fva_max_pw - fva_min_pw if all(np.isfinite([fva_min_pw, fva_max_pw])) else np.nan

		fva_pair_rows.append({
			"t_h":              t_fva,
			"mode":             mode_fva,
			"reaction":         rid,
			"label":            PAIR_LABELS.get(rid, rid),
			"fva_min_baseline": fva_min_b,
			"fva_max_baseline": fva_max_b,
			"fva_min_pairwise": fva_min_pw,
			"fva_max_pairwise": fva_max_pw,
			"range_baseline":   range_b,
			"range_pairwise":   range_pw,
			"delta_range":      range_pw - range_b if all(np.isfinite([range_b, range_pw])) else np.nan,
		})

ester_fva_pairwise_df = pd.DataFrame(fva_pair_rows)

print("\n[FVA baseline vs pairwise – 3 pares activos (fraction_of_optimum=0.95)]")
if not ester_fva_pairwise_df.empty:
	display(ester_fva_pairwise_df[[
		"t_h", "mode", "reaction", "label",
		"fva_min_baseline", "fva_max_baseline",
		"fva_min_pairwise", "fva_max_pairwise",
		"delta_range",
	]])
else:
	print("[INFO] FVA no ejecutada (sin snapshots en phase_times_fva o cache vacío).")

# ──── Verificación de signo: r_1865, r_1866, r_1589 en tiempos de fase ───────────────
print("\n" + "─" * 80)
print("Verificación de signo: alcoholes del par en tiempos clave")
print("  Esperado: flujo pairwise > 0  (exportación neta del metabolito al medio)")
print("─" * 80)

sign_rows = []

for t_sv in phase_times_fva:
	tk_row_sv = ester_pairwise_time_df[ester_pairwise_time_df["t_h"] == float(t_sv)]
	if tk_row_sv.empty:
		continue
	tr_sv   = tk_row_sv.iloc[0]
	mode_sv = tr_sv["mode"]

	for alc_rid_sv in ALCOHOL_SIGN_IDS:
		v_pw_sv = float(tr_sv.get(f"alcohol_{alc_rid_sv}_pairwise", np.nan))
		v_b_sv  = float(tr_sv.get(f"alcohol_{alc_rid_sv}_baseline",  np.nan))
		if np.isfinite(v_pw_sv) and abs(v_pw_sv) < ZERO_TOL: v_pw_sv = 0.0
		if np.isfinite(v_b_sv)  and abs(v_b_sv)  < ZERO_TOL: v_b_sv  = 0.0

		# También leer fva_max_pairwise del df de FVA si está disponible
		fva_max_pw_sv = np.nan
		if not ester_fva_pairwise_df.empty:
			_fva_row = ester_fva_pairwise_df[
				(ester_fva_pairwise_df["t_h"] == float(t_sv)) &
				(ester_fva_pairwise_df["reaction"] == alc_rid_sv)
			]
			if not _fva_row.empty:
				fva_max_pw_sv = float(_fva_row["fva_max_pairwise"].iloc[0])

		sign_ok  = (v_pw_sv > 0.0) if np.isfinite(v_pw_sv) else None
		sign_str = "OK  ✓" if sign_ok else ("cero/neg ✗" if sign_ok is False else "NaN  ?")
		label_sv = PAIR_LABELS.get(alc_rid_sv, alc_rid_sv)

		print(f"  t={t_sv:5.1f}h  [{mode_sv:12s}]  {alc_rid_sv}  ({label_sv:20s})"
		      f"  v_pw={v_pw_sv:+.4e}   v_b={v_b_sv:+.4e}"
		      f"   fva_max_pw={fva_max_pw_sv:+.4e}   {sign_str}")

		sign_rows.append({
			"t_h":          t_sv,
			"mode":         mode_sv,
			"alc_rid":      alc_rid_sv,
			"label":        label_sv,
			"v_baseline":   v_b_sv,
			"v_pairwise":   v_pw_sv,
			"fva_max_pw":   fva_max_pw_sv,
			"sign_ok":      sign_ok,
		})

ester_alcohol_sign_df = pd.DataFrame(sign_rows)
if not ester_alcohol_sign_df.empty:
	display(ester_alcohol_sign_df)

# ─────────────────────────────────────────────────────────────────────────────────────
# F. CONCENTRACIÓN MACROSCÓPICA ACUMULADA (pairwise)
# ─────────────────────────────────────────────────────────────────────────────────────

def _get_exchange_mw_g_per_mmol_pw(mdl, rid):
	rxn  = mdl.reactions.get_by_id(rid)
	mets = list(rxn.metabolites.keys())
	if len(mets) != 1:
		return np.nan
	met = mets[0]
	try:
		return float(met.formula_weight) / 1000.0
	except Exception:
		return np.nan

ester_mw_map_pw      = {rid: _get_exchange_mw_g_per_mmol_pw(model, rid) for rid in acetate_ester_targets}
ester_conc_baseline_pw = {rid: np.zeros(len(state_times_pw), dtype=float) for rid in acetate_ester_targets}
ester_conc_pairwise_pw = {rid: np.zeros(len(state_times_pw), dtype=float) for rid in acetate_ester_targets}

for k, t_k in enumerate(analysis_times_pw):
	row_state_e = states_df_phase.loc[states_df_phase["t_h"] == float(t_k)].iloc[0]
	Xk          = float(row_state_e["X_gDW_L"]) if "X_gDW_L" in row_state_e.index else 0.0

	idx_rows_e = ester_pairwise_time_df[ester_pairwise_time_df["t_h"] == float(t_k)].index
	if idx_rows_e.empty:
		continue
	idx_row_e = idx_rows_e[0]

	if k + 1 < len(state_times_pw):
		dt_k = float(state_times_pw[k + 1] - state_times_pw[k])
		for rid in acetate_ester_targets:
			mw      = ester_mw_map_pw.get(rid, np.nan)
			vb_col  = f"ester_{rid}_baseline"
			vpw_col = f"ester_{rid}_pairwise"
			vb_e  = float(ester_pairwise_time_df.loc[idx_row_e, vb_col])  if vb_col  in ester_pairwise_time_df.columns else np.nan
			vpw_e = float(ester_pairwise_time_df.loc[idx_row_e, vpw_col]) if vpw_col in ester_pairwise_time_df.columns else np.nan
			vb_e  = max(0.0, vb_e)  if np.isfinite(vb_e)  else 0.0
			vpw_e = max(0.0, vpw_e) if np.isfinite(vpw_e) else 0.0
			ester_conc_baseline_pw[rid][k + 1] = ester_conc_baseline_pw[rid][k]
			ester_conc_pairwise_pw[rid][k + 1] = ester_conc_pairwise_pw[rid][k]
			if np.isfinite(mw):
				ester_conc_baseline_pw[rid][k + 1] += 1000.0 * mw * vb_e  * Xk * dt_k
				ester_conc_pairwise_pw[rid][k + 1] += 1000.0 * mw * vpw_e * Xk * dt_k

ester_pairwise_conc_df = pd.DataFrame({"t_h": state_times_pw})
for rid, label in acetate_ester_targets.items():
	sfx = re.sub(r"[^0-9a-zA-Z]+", "_", str(label)).lower()
	ester_pairwise_conc_df[f"{sfx}_baseline_mg_L"] = ester_conc_baseline_pw[rid]
	ester_pairwise_conc_df[f"{sfx}_pairwise_mg_L"] = ester_conc_pairwise_pw[rid]
	ester_pairwise_conc_df[f"{sfx}_delta_mg_L"]    = ester_conc_pairwise_pw[rid] - ester_conc_baseline_pw[rid]

baseline_conc_cols_pw = [c for c in ester_pairwise_conc_df.columns if c.endswith("_baseline_mg_L")]
pairwise_conc_cols_pw = [c for c in ester_pairwise_conc_df.columns if c.endswith("_pairwise_mg_L")]
ester_pairwise_conc_df["total_acetate_esters_baseline_mg_L"] = ester_pairwise_conc_df[baseline_conc_cols_pw].sum(axis=1)
ester_pairwise_conc_df["total_acetate_esters_pairwise_mg_L"] = ester_pairwise_conc_df[pairwise_conc_cols_pw].sum(axis=1)
ester_pairwise_conc_df["total_acetate_esters_delta_mg_L"]    = (
	ester_pairwise_conc_df["total_acetate_esters_pairwise_mg_L"]
	- ester_pairwise_conc_df["total_acetate_esters_baseline_mg_L"]
)

print("\n[Concentración acumulada: baseline vs pairwise]")
display(ester_pairwise_conc_df[[
	"t_h",
	"total_acetate_esters_baseline_mg_L",
	"total_acetate_esters_pairwise_mg_L",
	"total_acetate_esters_delta_mg_L"
]])

# ─────────────────────────────────────────────────────────────────────────────────────
# G. GRÁFICOS
# ─────────────────────────────────────────────────────────────────────────────────────

import matplotlib.pyplot as plt

fig, ax = plt.subplots(2, 2, figsize=(14, 8), sharex=True)

ax[0, 0].plot(ester_pairwise_conc_df["t_h"], ester_pairwise_conc_df["total_acetate_esters_baseline_mg_L"],
              lw=2.5, label="baseline", color="tab:blue")
ax[0, 0].plot(ester_pairwise_conc_df["t_h"], ester_pairwise_conc_df["total_acetate_esters_pairwise_mg_L"],
              lw=2.5, ls="--", label="pairwise", color="tab:orange")
ax[0, 0].set_title("Total acetate esters (mg/L) – Pairwise")
ax[0, 0].set_ylabel("mg/L"); ax[0, 0].grid(alpha=0.3); ax[0, 0].legend()

active_esters_for_plot = list(active_pairs["ester_rid"].unique())[:4]
colors_esters = plt.cm.tab10(np.linspace(0, 1, max(len(active_esters_for_plot), 1)))

for i, ester_rid_p in enumerate(active_esters_for_plot):
	label_e = acetate_ester_targets.get(ester_rid_p, ester_rid_p)
	sfx     = re.sub(r"[^0-9a-zA-Z]+", "_", str(label_e)).lower()
	col_b   = f"{sfx}_baseline_mg_L"
	col_pw  = f"{sfx}_pairwise_mg_L"
	if col_b in ester_pairwise_conc_df.columns and col_pw in ester_pairwise_conc_df.columns:
		ax[0, 1].plot(ester_pairwise_conc_df["t_h"], ester_pairwise_conc_df[col_b],
		              lw=1.8, color=colors_esters[i], alpha=0.6, label=f"{label_e} (bsl)")
		ax[0, 1].plot(ester_pairwise_conc_df["t_h"], ester_pairwise_conc_df[col_pw],
		              lw=1.8, ls="--", color=colors_esters[i], alpha=0.9, label=f"{label_e} (pair)")
ax[0, 1].set_title("Ésters activos en pairwise (mg/L)")
ax[0, 1].set_ylabel("mg/L"); ax[0, 1].grid(alpha=0.3); ax[0, 1].legend(fontsize=8, loc="best")

ax[1, 0].plot(ester_pairwise_time_df["t_h"], ester_pairwise_time_df["mu_h_inv_baseline"],
              lw=2, label="mu baseline")
ax[1, 0].plot(ester_pairwise_time_df["t_h"], ester_pairwise_time_df["mu_h_inv_pairwise"],
              lw=2, ls="--", label="mu pairwise")
ax[1, 0].set_title("μ (h⁻¹) – Baseline vs Pairwise")
ax[1, 0].set_xlabel("t (h)"); ax[1, 0].set_ylabel("μ (h⁻¹)")
ax[1, 0].grid(alpha=0.3); ax[1, 0].legend()

active_pairs_list = active_pairs.to_dict("records")
if active_pairs_list:
	binding_counts = [
		sum(1 for pair in active_pairs_list
		    if row_df.get(f"binding_{pair['ester_rid']}__{pair['alcohol_rid']}", False))
		for _, row_df in ester_pairwise_time_df.iterrows()
	]
	ax[1, 1].bar(ester_pairwise_time_df["t_h"], binding_counts, width=0.8, alpha=0.7, color="tab:red")
	ax[1, 1].axhline(len(active_pairs), ls="--", lw=1, color="gray", label="max activas")
	ax[1, 1].set_title("# restricciones pairwise activas (binding)")
	ax[1, 1].set_xlabel("t (h)"); ax[1, 1].set_ylabel("# binding")
	ax[1, 1].grid(alpha=0.3); ax[1, 1].legend()

plt.suptitle("Comparativa pairwise: Reforma de acoplamiento ester↔alcohol", fontsize=12, y=1.01)
plt.tight_layout()
plt.show()

# ─────────────────────────────────────────────────────────────────────────────────────
# H. SENSIBILIDAD EN PHI_PAIR
# ─────────────────────────────────────────────────────────────────────────────────────

ester_pairwise_sensitivity_df = pd.DataFrame()

if RUN_PAIRWISE_SENSITIVITY:
	sens_rows_pw = []
	phi_grid_pw  = sorted(set(float(x) for x in PAIRWISE_PHI_GRID if float(x) >= 0.0))
	phase_times_pw = set(float(x) for x in phase_meta_df["t_h"].tolist()) if not phase_meta_df.empty else set()

	for t_phase_pk in phase_times_pw:
		for phi_pair_val in phi_grid_pw:
			pairwise_df_phi = ester_alcohol_pairwise_df.copy()
			pairwise_df_phi.loc[pairwise_df_phi["pair_active"], "phi_pair"] = float(phi_pair_val)

			res_phi    = _build_phase_snapshot_model_pairwise(
				t_phase_pk,
				enable_pairwise=(phi_pair_val > 0.0),
				pairwise_df=pairwise_df_phi,
			)
			sol_phi    = res_phi["solution_pairwise"]
			status_phi = getattr(sol_phi, "status", "unknown")

			sens_rows_pw.append({
				"t_h":             float(t_phase_pk),
				"phi_pair":        phi_pair_val,
				"status":          status_phi,
				"objective_id":    res_phi["objective_id"],
				"objective_value": float(sol_phi.fluxes.get(res_phi["objective_id"], np.nan)) if status_phi == "optimal" else np.nan,
				"mu_h_inv":        float(sol_phi.fluxes.get(OBJ_ID, np.nan)) if status_phi == "optimal" else np.nan,
				"n_binding":       sum(1 for v in res_phi["constraints_binding"].values() if v["binding"]),
			})

	ester_pairwise_sensitivity_df = pd.DataFrame(sens_rows_pw)
	if not ester_pairwise_sensitivity_df.empty:
		print("\n[Sensibilidad pairwise en phi_pair]")
		display(ester_pairwise_sensitivity_df)

# ─────────────────────────────────────────────────────────────────────────────────────
# I. RESUMEN FINAL Y GUARDADO
# ─────────────────────────────────────────────────────────────────────────────────────

summary_rows_pw = []
for rid, label in acetate_ester_targets.items():
	sfx    = re.sub(r"[^0-9a-zA-Z]+", "_", str(label)).lower()
	col_b  = f"{sfx}_baseline_mg_L"
	col_pw = f"{sfx}_pairwise_mg_L"
	final_b  = float(ester_pairwise_conc_df[col_b].iloc[-1])  if col_b  in ester_pairwise_conc_df.columns else 0.0
	final_pw = float(ester_pairwise_conc_df[col_pw].iloc[-1]) if col_pw in ester_pairwise_conc_df.columns else 0.0
	is_in_pair = len(active_pairs[active_pairs["ester_rid"] == rid]) > 0
	summary_rows_pw.append({
		"rid": rid, "label": label,
		"in_pairwise_active":  is_in_pair,
		"final_baseline_mg_L": final_b,
		"final_pairwise_mg_L": final_pw,
		"delta_mg_L":          final_pw - final_b,
	})

ester_pairwise_summary_df = pd.DataFrame(summary_rows_pw).sort_values(
	"delta_mg_L", ascending=False, kind="stable"
)
print("\n[Resumen final: comparativa pairwise vs baseline]")
display(ester_pairwise_summary_df)

if SAVE_PAIRWISE_CSV:
	pairwise_map_file     = OUT_DIR / "ester_pairwise_mapping_cell25.csv"
	pairwise_time_file    = OUT_DIR / "ester_pairwise_time_comparison_cell25.csv"
	pairwise_conc_file    = OUT_DIR / "ester_pairwise_concentration_comparison_cell25.csv"
	pairwise_summary_file = OUT_DIR / "ester_pairwise_summary_cell25.csv"
	pairwise_diag_file    = OUT_DIR / "ester_pairwise_diagnostic_cell25.csv"
	pairwise_fva_file     = OUT_DIR / "ester_pairwise_fva_cell25.csv"
	pairwise_sign_file    = OUT_DIR / "ester_pairwise_sign_check_cell25.csv"

	ester_alcohol_pairwise_df.to_csv(pairwise_map_file,                  index=False)
	ester_pairwise_time_df.to_csv(pairwise_time_file,                    index=False)
	ester_pairwise_conc_df.to_csv(pairwise_conc_file,                    index=False)
	ester_pairwise_summary_df.to_csv(pairwise_summary_file,              index=False)
	ester_pairwise_pair_diagnostics_df.to_csv(pairwise_diag_file,        index=False)
	if not ester_fva_pairwise_df.empty:
		ester_fva_pairwise_df.to_csv(pairwise_fva_file,                  index=False)
	if not ester_alcohol_sign_df.empty:
		ester_alcohol_sign_df.to_csv(pairwise_sign_file,                 index=False)

	print("\nArchivos guardados:")
	for f_path in [pairwise_map_file, pairwise_time_file, pairwise_conc_file,
	               pairwise_summary_file, pairwise_diag_file, pairwise_fva_file,
	               pairwise_sign_file]:
		print(f"  {f_path}")

	if not ester_pairwise_sensitivity_df.empty:
		pairwise_sens_file = OUT_DIR / "ester_pairwise_sensitivity_cell25.csv"
		ester_pairwise_sensitivity_df.to_csv(pairwise_sens_file, index=False)
		print(f"  {pairwise_sens_file}")

print("\n[Variables disponibles]")
for vname in [
	"ester_alcohol_pairwise_df",          "ester_pairwise_time_df",
	"ester_pairwise_conc_df",             "ester_pairwise_summary_df",
	"ester_pairwise_pair_diagnostics_df", "ester_fva_pairwise_df",
	"ester_alcohol_sign_df",              "ester_pairwise_sensitivity_df",
	"snapshot_cache_pairwise",
]:
	print(f"  {vname}")

print("\n" + "=" * 100)
print("OK: Reforma pairwise completada (snapshot-level sobre trayectoria baseline).")
print("IDs manuales corregidos: r_1862->r_1865 | r_1867->r_1866 | r_2000->r_1589")
print("ethyl acetate (r_1765) excluido del acoplamiento pairwise principal.")
print("=" * 100)
