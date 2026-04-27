# ======================================================================================
# Acetate esters: REINTRODUCCIÓN SUAVE DE ETHYL ACETATE (soft coupling independiente)
# Análisis snapshot-level sobre trayectoria baseline – INCREMENTAL sobre celda anterior
# ======================================================================================

print("\n" + "=" * 100)
print("Ethyl acetate reintroducido con soft coupling independiente (snapshot-level)")
print("  [Incremento sobre pairwise: mantiene r_1867/r_1862/r_2000; agrega r_1765 por separado]")
print("=" * 100)

# ── Verificar dependencias ─────────────────────────────────────────────────────────
_required_sym_ea = [
    "model", "OUT_DIR",
    "states_df_phase", "flux_df_phase", "phase_meta_df",
    "_build_phase_snapshot_model",
    "acetate_ester_targets", "available_aroma_targets",
    "ester_alcohol_pairwise_df", "active_pairs",
    "OBJ_ID", "ATPM_ID", "ETH_ID",
    "ZERO_TOL", "ESTER_BIND_TOL",
    "PAIRWISE_PHI_STATIC",
    "_add_pairwise_ester_coupling_constraints",
    "_build_phase_snapshot_model_pairwise",
    "ester_pairwise_time_df",
    "ester_pairwise_conc_df",
    "ester_mw_map_pw",
]
_missing_ea = [s for s in _required_sym_ea if s not in globals()]
if _missing_ea:
    raise RuntimeError(f"Faltan símbolos requeridos (ejecuta celda anterior primero): {_missing_ea}")

# ── Configuración: soft coupling para ethyl acetate ───────────────────────────────
INCLUDE_ETHYL_ACETATE_SOFT_COUPLING = True
ETHYL_ACETATE_MODE                  = "soft_ethanol_coupling"  # "soft_ethanol_coupling" | "off"

ETHYL_ACETATE_RID_EA = "r_1765"   # ethyl acetate exchange
ETHANOL_RID_EA       = ETH_ID     # r_1761, ethanol exchange

PHI_ETHYL_ACETATE_STATIC = 0.0001  # ← PRIMER PARÁMETRO A TOCAR
PHI_ETHYL_ACETATE_GRID   = [0.000, 0.002, 0.005, 0.010, 0.020]

SAVE_EA_CSV = True

print(f"INCLUDE_ETHYL_ACETATE_SOFT_COUPLING = {INCLUDE_ETHYL_ACETATE_SOFT_COUPLING}")
print(f"ETHYL_ACETATE_MODE                  = {ETHYL_ACETATE_MODE}")
print(f"PHI_ETHYL_ACETATE_STATIC            = {PHI_ETHYL_ACETATE_STATIC}")
print(f"PHI_ETHYL_ACETATE_GRID              = {PHI_ETHYL_ACETATE_GRID}")
print(f"Ethyl acetate RID : {ETHYL_ACETATE_RID_EA}")
print(f"Ethanol RID       : {ETHANOL_RID_EA}")

# ── Verificar presencia en modelo ────────────────────────────────────────────────
_rxn_ids_all_ea = {r.id for r in model.reactions}
if ETHYL_ACETATE_RID_EA not in _rxn_ids_all_ea:
    raise RuntimeError(f"{ETHYL_ACETATE_RID_EA} no encontrado en modelo.")
if ETHANOL_RID_EA not in _rxn_ids_all_ea:
    raise RuntimeError(f"{ETHANOL_RID_EA} no encontrado en modelo.")

print(f"\n[OK] {ETHYL_ACETATE_RID_EA}: {model.reactions.get_by_id(ETHYL_ACETATE_RID_EA).name}")
print(f"[OK] {ETHANOL_RID_EA}: {model.reactions.get_by_id(ETHANOL_RID_EA).name}")

# ─────────────────────────────────────────────────────────────────────────────────
# A. WRAPPER: snapshot + pairwise + soft ethyl acetate coupling
# ─────────────────────────────────────────────────────────────────────────────────

def _add_ethyl_acetate_soft_constraint(mtmp, phi_ea):
    """
    Agrega restricción suave independiente:
        v_ethyl_acetate >= phi_ea * v_ethanol
      ⟺ v_ethyl_acetate - phi_ea * v_ethanol >= 0

    Retorna (constraint_object | None, added: bool)
    """
    if (not INCLUDE_ETHYL_ACETATE_SOFT_COUPLING
            or ETHYL_ACETATE_MODE == "off"
            or phi_ea <= 0.0):
        return None, False

    rb = {r.id: r for r in mtmp.reactions}
    if ETHYL_ACETATE_RID_EA not in rb or ETHANOL_RID_EA not in rb:
        return None, False

    ea_expr  = rb[ETHYL_ACETATE_RID_EA].flux_expression
    eth_expr = rb[ETHANOL_RID_EA].flux_expression

    cname = f"ethyl_acetate_soft_phi{str(phi_ea).replace('.', '_')}"
    cons  = mtmp.problem.Constraint(
        ea_expr - float(phi_ea) * eth_expr,
        lb=0.0,
        name=cname,
    )
    mtmp.add_cons_vars([cons])
    return cons, True


def _build_snapshot_pairwise_plus_ea(t_phase, *, phi_ea=None, pairwise_df=None):
    """
    Snapshot sobre trayectoria baseline con:
      1) restricciones pairwise actuales (r_1867/r_1862/r_2000)
      2) soft coupling independiente para ethyl acetate

    Retorna dict extendido con métricas de ethyl acetate.
    """
    phi_ea_val = float(PHI_ETHYL_ACETATE_STATIC if phi_ea is None else phi_ea)
    pw_df      = ester_alcohol_pairwise_df if pairwise_df is None else pairwise_df

    # Base snapshot (sin restricciones pairwise)
    snap = _build_phase_snapshot_model(float(t_phase))
    m_base  = snap["model"]
    sol_base = snap["solution"]

    # Modelo con pairwise + ethyl acetate soft
    m_full = m_base.copy()

    # 1) Pairwise principal
    _, active_pair_map, phi_map = _add_pairwise_ester_coupling_constraints(m_full, pw_df)

    # 2) Soft coupling para ethyl acetate
    ea_cons, ea_added = _add_ethyl_acetate_soft_constraint(m_full, phi_ea_val)

    try:
        sol_full = m_full.optimize()
    except Exception as exc:
        print(f"[WARN] t={t_phase}: optimize falló: {exc}")
        sol_full = sol_base

    status_full = getattr(sol_full, "status", "unknown")

    def _safe(sol, rid):
        if getattr(sol, "status", "") != "optimal":
            return np.nan
        v = float(sol.fluxes.get(rid, 0.0))
        return 0.0 if abs(v) < ZERO_TOL else v

    # Binding pairwise
    pair_binding = {}
    if status_full == "optimal":
        for (e_rid, a_rid) in active_pair_map:
            phi_p = phi_map.get((e_rid, a_rid), 0.0)
            v_e   = _safe(sol_full, e_rid)
            v_a   = _safe(sol_full, a_rid)
            lhs   = v_e - phi_p * v_a if all(np.isfinite([v_e, v_a])) else np.nan
            if np.isfinite(lhs) and abs(lhs) < ZERO_TOL:
                lhs = 0.0
            pair_binding[(e_rid, a_rid)] = {
                "phi_pair": phi_p, "v_ester": v_e, "v_alcohol": v_a,
                "lhs_margin": lhs, "binding": lhs <= ESTER_BIND_TOL if np.isfinite(lhs) else False,
            }

    # Ethyl acetate binding
    v_ea_b   = _safe(sol_base, ETHYL_ACETATE_RID_EA)
    v_ea_f   = _safe(sol_full, ETHYL_ACETATE_RID_EA)
    v_eth_b  = _safe(sol_base, ETHANOL_RID_EA)
    v_eth_f  = _safe(sol_full, ETHANOL_RID_EA)

    lhs_ea   = v_ea_f - phi_ea_val * v_eth_f if all(np.isfinite([v_ea_f, v_eth_f])) else np.nan
    if np.isfinite(lhs_ea) and abs(lhs_ea) < ZERO_TOL:
        lhs_ea = 0.0
    binding_ea = (lhs_ea <= ESTER_BIND_TOL) if np.isfinite(lhs_ea) else False

    return {
        "t_h":              float(t_phase),
        "mode":             snap["mode"],
        "objective_id":     snap["objective_id"],
        "state_row":        snap["state_row"],
        "model_baseline":   m_base,
        "model_full":       m_full,
        "solution_baseline": sol_base,
        "solution_full":     sol_full,
        "phi_ea":            phi_ea_val,
        "ea_constraint_added": ea_added,
        "pair_binding":      pair_binding,
        # ethyl acetate
        "v_ethanol_baseline":       v_eth_b,
        "v_ethanol_full":           v_eth_f,
        "v_eth_acetate_baseline":   v_ea_b,
        "v_eth_acetate_full":       v_ea_f,
        "lhs_margin_eth":           lhs_ea,
        "binding_eth":              binding_ea,
    }


# ─────────────────────────────────────────────────────────────────────────────────
# B. COMPARACIÓN TEMPORAL baseline vs (pairwise + soft EA)
# ─────────────────────────────────────────────────────────────────────────────────

analysis_times_ea = flux_df_phase["t_h"].to_numpy(dtype=float)
state_times_ea    = states_df_phase["t_h"].to_numpy(dtype=float)

ea_comparison_rows = []
snapshot_cache_ea  = {}

for k, t_k in enumerate(analysis_times_ea):
    res_ea = _build_snapshot_pairwise_plus_ea(t_k)
    snapshot_cache_ea[float(t_k)] = {
        "model":    res_ea["model_full"],
        "solution": res_ea["solution_full"],
    }

    sol_b  = res_ea["solution_baseline"]
    sol_f  = res_ea["solution_full"]
    obj_id = res_ea["objective_id"]

    def _sf(sol, rid):
        if getattr(sol, "status", "") != "optimal":
            return np.nan
        v = float(sol.fluxes.get(rid, 0.0))
        return 0.0 if abs(v) < ZERO_TOL else v

    row = {
        "t_h":                         float(t_k),
        "mode":                        res_ea["mode"],
        "objective_id":                obj_id,
        "objective_value_baseline":    _sf(sol_b, obj_id),
        "objective_value_full":        _sf(sol_f, obj_id),
        "delta_objective":             (_sf(sol_f, obj_id) - _sf(sol_b, obj_id))
                                       if all(np.isfinite([_sf(sol_f, obj_id), _sf(sol_b, obj_id)])) else np.nan,
        "mu_baseline":                 _sf(sol_b, OBJ_ID),
        "mu_full":                     _sf(sol_f, OBJ_ID),
        "v_ethanol_baseline":          res_ea["v_ethanol_baseline"],
        "v_ethanol_full":              res_ea["v_ethanol_full"],
        "v_eth_acetate_baseline":      res_ea["v_eth_acetate_baseline"],
        "v_eth_acetate_full":          res_ea["v_eth_acetate_full"],
        "phi_ea":                      res_ea["phi_ea"],
        "lhs_margin_eth":              res_ea["lhs_margin_eth"],
        "binding_eth":                 res_ea["binding_eth"],
        "ea_constraint_added":         res_ea["ea_constraint_added"],
        "n_pairwise_binding":          sum(1 for v in res_ea["pair_binding"].values() if v["binding"]),
    }

    # Flujos individuales de todos los acetate esters
    for e_rid, e_lbl in acetate_ester_targets.items():
        row[f"ester_{e_rid}_baseline"] = _sf(sol_b, e_rid)
        row[f"ester_{e_rid}_full"]     = _sf(sol_f, e_rid)

    ea_comparison_rows.append(row)

ea_comparison_df = pd.DataFrame(ea_comparison_rows)

print("\n[Diagnóstico ethyl acetate – primeros y últimos 5 pasos]")
_diag_cols_ea = [
    "t_h", "mode",
    "v_ethanol_baseline", "v_ethanol_full",
    "v_eth_acetate_baseline", "v_eth_acetate_full",
    "phi_ea", "lhs_margin_eth", "binding_eth",
]
display(pd.concat([ea_comparison_df[_diag_cols_ea].head(5),
                   ea_comparison_df[_diag_cols_ea].tail(5)]))

# ─────────────────────────────────────────────────────────────────────────────────
# C. CONCENTRACIÓN MACROSCÓPICA ACUMULADA
# ─────────────────────────────────────────────────────────────────────────────────

ester_conc_ea_baseline = {rid: np.zeros(len(state_times_ea), dtype=float) for rid in acetate_ester_targets}
ester_conc_ea_full     = {rid: np.zeros(len(state_times_ea), dtype=float) for rid in acetate_ester_targets}

for k, t_k in enumerate(analysis_times_ea):
    if k + 1 >= len(state_times_ea):
        continue

    row_state = states_df_phase.loc[states_df_phase["t_h"] == float(t_k)].iloc[0]
    Xk        = float(row_state["X_gDW_L"]) if "X_gDW_L" in row_state.index else 0.0
    dt_k      = float(state_times_ea[k + 1] - state_times_ea[k])

    idx_rows  = ea_comparison_df[ea_comparison_df["t_h"] == float(t_k)].index
    if idx_rows.empty:
        continue
    idx_row = idx_rows[0]

    for rid in acetate_ester_targets:
        mw      = ester_mw_map_pw.get(rid, np.nan)
        if not np.isfinite(mw):
            continue
        col_b   = f"ester_{rid}_baseline"
        col_f   = f"ester_{rid}_full"
        vb_k    = max(0.0, float(ea_comparison_df.loc[idx_row, col_b])) if col_b in ea_comparison_df.columns and np.isfinite(ea_comparison_df.loc[idx_row, col_b]) else 0.0
        vf_k    = max(0.0, float(ea_comparison_df.loc[idx_row, col_f])) if col_f in ea_comparison_df.columns and np.isfinite(ea_comparison_df.loc[idx_row, col_f]) else 0.0

        ester_conc_ea_baseline[rid][k + 1] = ester_conc_ea_baseline[rid][k] + 1000.0 * mw * vb_k * Xk * dt_k
        ester_conc_ea_full[rid][k + 1]     = ester_conc_ea_full[rid][k]     + 1000.0 * mw * vf_k * Xk * dt_k

ea_conc_df = pd.DataFrame({"t_h": state_times_ea})
for rid, lbl in acetate_ester_targets.items():
    sfx = re.sub(r"[^0-9a-zA-Z]+", "_", str(lbl)).lower()
    ea_conc_df[f"{sfx}_baseline_mg_L"] = ester_conc_ea_baseline[rid]
    ea_conc_df[f"{sfx}_full_mg_L"]     = ester_conc_ea_full[rid]

b_cols = [c for c in ea_conc_df.columns if c.endswith("_baseline_mg_L")]
f_cols = [c for c in ea_conc_df.columns if c.endswith("_full_mg_L")]
ea_conc_df["total_baseline_mg_L"] = ea_conc_df[b_cols].sum(axis=1)
ea_conc_df["total_full_mg_L"]     = ea_conc_df[f_cols].sum(axis=1)
ea_conc_df["total_delta_mg_L"]    = ea_conc_df["total_full_mg_L"] - ea_conc_df["total_baseline_mg_L"]

# ─────────────────────────────────────────────────────────────────────────────────
# D. RESUMEN DISTRIBUCIÓN DEL POOL DE ÉSTERES (final t=72 h)
# ─────────────────────────────────────────────────────────────────────────────────

pool_rows = []
total_final_b  = float(ea_conc_df["total_baseline_mg_L"].iloc[-1])
total_final_f  = float(ea_conc_df["total_full_mg_L"].iloc[-1])

for rid, lbl in acetate_ester_targets.items():
    sfx   = re.sub(r"[^0-9a-zA-Z]+", "_", str(lbl)).lower()
    col_b = f"{sfx}_baseline_mg_L"
    col_f = f"{sfx}_full_mg_L"
    fin_b = float(ea_conc_df[col_b].iloc[-1]) if col_b in ea_conc_df.columns else 0.0
    fin_f = float(ea_conc_df[col_f].iloc[-1]) if col_f in ea_conc_df.columns else 0.0
    pool_rows.append({
        "rid":                   rid,
        "label":                 lbl,
        "final_baseline_mg_L":  fin_b,
        "final_full_mg_L":      fin_f,
        "delta_mg_L":           fin_f - fin_b,
        "frac_baseline_%":      100.0 * fin_b / total_final_b if total_final_b > 0 else 0.0,
        "frac_full_%":          100.0 * fin_f / total_final_f if total_final_f > 0 else 0.0,
    })

pool_rows.append({
    "rid": "TOTAL", "label": "Total acetate esters",
    "final_baseline_mg_L": total_final_b,
    "final_full_mg_L":     total_final_f,
    "delta_mg_L":          total_final_f - total_final_b,
    "frac_baseline_%":     100.0,
    "frac_full_%":         100.0,
})

ea_pool_df = pd.DataFrame(pool_rows)
print("\n[Distribución del pool de acetate esters (t=72 h)]")
display(ea_pool_df)

# ─────────────────────────────────────────────────────────────────────────────────
# E. SENSIBILIDAD EN PHI_ETHYL_ACETATE
# ─────────────────────────────────────────────────────────────────────────────────

phi_grid_ea    = sorted(set(float(x) for x in PHI_ETHYL_ACETATE_GRID))
phase_times_ea_set = sorted(set(float(x) for x in phase_meta_df["t_h"].tolist())) if not phase_meta_df.empty else []

sens_ea_rows = []

for t_phase_ea in phase_times_ea_set:
    for phi_val_ea in phi_grid_ea:
        res_phi_ea = _build_snapshot_pairwise_plus_ea(t_phase_ea, phi_ea=phi_val_ea)
        sol_phi_ea = res_phi_ea["solution_full"]
        status_ea  = getattr(sol_phi_ea, "status", "unknown")

        # Concentración final extrapolada: integrar desde t_phase_ea a t=72 con flujo constante
        Xk_ea  = float(res_phi_ea["state_row"].get("X_gDW_L", 1.0)) if hasattr(res_phi_ea["state_row"], "get") else float(res_phi_ea["state_row"]["X_gDW_L"]) if "X_gDW_L" in res_phi_ea["state_row"].index else 1.0
        dt_rem = max(0.0, 72.0 - float(t_phase_ea))

        def _sf2(sol, rid):
            if getattr(sol, "status", "") != "optimal":
                return np.nan
            v = float(sol.fluxes.get(rid, 0.0))
            return 0.0 if abs(v) < ZERO_TOL else v

        v_ea_phi = max(0.0, _sf2(sol_phi_ea, ETHYL_ACETATE_RID_EA)) if np.isfinite(_sf2(sol_phi_ea, ETHYL_ACETATE_RID_EA)) else 0.0
        mw_ea    = ester_mw_map_pw.get(ETHYL_ACETATE_RID_EA, np.nan)
        ea_mg_L_extrap = 1000.0 * mw_ea * v_ea_phi * Xk_ea * dt_rem if np.isfinite(mw_ea) else np.nan

        # Total ésteres extrapolado
        total_extrap = 0.0
        for rid_s in acetate_ester_targets:
            mw_s  = ester_mw_map_pw.get(rid_s, np.nan)
            v_s   = max(0.0, _sf2(sol_phi_ea, rid_s)) if np.isfinite(_sf2(sol_phi_ea, rid_s)) else 0.0
            if np.isfinite(mw_s):
                total_extrap += 1000.0 * mw_s * v_s * Xk_ea * dt_rem

        frac_ea = 100.0 * ea_mg_L_extrap / total_extrap if (total_extrap > 0 and np.isfinite(ea_mg_L_extrap)) else np.nan

        sens_ea_rows.append({
            "t_h":                     float(t_phase_ea),
            "mode":                    res_phi_ea["mode"],
            "phi_ea":                  phi_val_ea,
            "status":                  status_ea,
            "objective_value":         _sf2(sol_phi_ea, res_phi_ea["objective_id"]),
            "mu_h_inv":                _sf2(sol_phi_ea, OBJ_ID),
            "v_eth_acetate_mmol_gDW_h": v_ea_phi,
            "v_ethanol_mmol_gDW_h":    _sf2(sol_phi_ea, ETHANOL_RID_EA),
            "ea_mg_L_extrap":          ea_mg_L_extrap,
            "total_ester_mg_L_extrap": total_extrap if total_extrap > 0 else np.nan,
            "frac_ea_pct":             frac_ea,
            "lhs_margin_eth":          res_phi_ea["lhs_margin_eth"],
            "binding_eth":             res_phi_ea["binding_eth"],
            "ea_constraint_added":     res_phi_ea["ea_constraint_added"],
        })

ea_sensitivity_df = pd.DataFrame(sens_ea_rows)
print("\n[Sensibilidad phi_eth por fase]")
display(ea_sensitivity_df[[
    "t_h", "mode", "phi_ea", "status",
    "objective_value", "mu_h_inv",
    "v_eth_acetate_mmol_gDW_h", "ea_mg_L_extrap",
    "total_ester_mg_L_extrap", "frac_ea_pct",
]])

# ─────────────────────────────────────────────────────────────────────────────────
# F. GRÁFICOS
# ─────────────────────────────────────────────────────────────────────────────────

fig, ax = plt.subplots(2, 3, figsize=(17, 9), dpi=120)

# F1. Total acetate esters
ax[0, 0].plot(ea_conc_df["t_h"], ea_conc_df["total_baseline_mg_L"],
              lw=2.5, label="baseline", color="tab:blue")
ax[0, 0].plot(ea_conc_df["t_h"], ea_conc_df["total_full_mg_L"],
              lw=2.5, ls="--", label=f"pairwise + soft EA (φ={PHI_ETHYL_ACETATE_STATIC})", color="tab:orange")
ax[0, 0].set_title("Total acetate esters (mg/L)")
ax[0, 0].set_ylabel("mg/L")
ax[0, 0].grid(alpha=0.3)
ax[0, 0].legend(fontsize=8)

# F2. Ethyl acetate individual
ea_sfx_ea = re.sub(r"[^0-9a-zA-Z]+", "_", "ethyl acetate exchange").lower()
if f"{ea_sfx_ea}_baseline_mg_L" in ea_conc_df.columns:
    ax[0, 1].plot(ea_conc_df["t_h"], ea_conc_df[f"{ea_sfx_ea}_baseline_mg_L"],
                  lw=2.5, label="baseline", color="tab:blue")
    ax[0, 1].plot(ea_conc_df["t_h"], ea_conc_df[f"{ea_sfx_ea}_full_mg_L"],
                  lw=2.5, ls="--", label="soft EA", color="tab:orange")
ax[0, 1].set_title("Ethyl acetate (mg/L)")
ax[0, 1].grid(alpha=0.3)
ax[0, 1].legend(fontsize=8)
ax[0, 1].set_ylabel("mg/L")

# F3. Distribución del pool (full) – barras
pool_plot = ea_pool_df[ea_pool_df["rid"] != "TOTAL"].copy()
valid_pool = pool_plot[pool_plot["final_full_mg_L"] > 1e-6]
if not valid_pool.empty:
    ax[0, 2].barh(valid_pool["label"], valid_pool["final_full_mg_L"], color="tab:green", alpha=0.8)
    ax[0, 2].set_title("Pool ésteres final (pairwise+soft EA, mg/L)")
    ax[0, 2].set_xlabel("mg/L")
    ax[0, 2].grid(alpha=0.3, axis="x")
else:
    ax[0, 2].text(0.5, 0.5, "Sin flujo detectado", ha="center", va="center")
    ax[0, 2].set_title("Pool ésteres final")

# F4. Diagnóstico ethyl acetate: v_eth_acetate vs t
ax[1, 0].plot(ea_comparison_df["t_h"], ea_comparison_df["v_eth_acetate_baseline"],
              lw=2, label="v_ea baseline", color="tab:blue")
ax[1, 0].plot(ea_comparison_df["t_h"], ea_comparison_df["v_eth_acetate_full"],
              lw=2, ls="--", label="v_ea full", color="tab:orange")
ax[1, 0].set_title("Flujo ethyl acetate (mmol/gDW/h)")
ax[1, 0].set_xlabel("t (h)")
ax[1, 0].set_ylabel("mmol/gDW/h")
ax[1, 0].grid(alpha=0.3)
ax[1, 0].legend(fontsize=8)

# F5. Fracción relativa de ethyl acetate (full) vs t
_ea_full_conc    = ea_conc_df[f"{ea_sfx_ea}_full_mg_L"].to_numpy(dtype=float) if f"{ea_sfx_ea}_full_mg_L" in ea_conc_df.columns else np.zeros(len(state_times_ea))
_total_full_conc = ea_conc_df["total_full_mg_L"].to_numpy(dtype=float)
_frac_ea_t       = np.where(_total_full_conc > 0, 100.0 * _ea_full_conc / _total_full_conc, 0.0)
ax[1, 1].plot(ea_conc_df["t_h"], _frac_ea_t, lw=2.5, color="tab:red")
ax[1, 1].set_title("Fracción ethyl acetate / total ésteres (%)")
ax[1, 1].set_xlabel("t (h)")
ax[1, 1].set_ylabel("%")
ax[1, 1].set_ylim(0, 105)
ax[1, 1].grid(alpha=0.3)

# F6. Sensibilidad phi_ea: frac_ea_pct vs phi en fase estacionaria
if not ea_sensitivity_df.empty:
    stat_sens = ea_sensitivity_df[ea_sensitivity_df["mode"] == "TURNOVER_ATPM"].copy()
    if stat_sens.empty:
        stat_sens = ea_sensitivity_df.copy()
    for t_ph, grp in stat_sens.groupby("t_h"):
        ax[1, 2].plot(grp["phi_ea"], grp["frac_ea_pct"], marker="o", lw=2, label=f"t={t_ph:.0f}h")
    ax[1, 2].set_title("Fracción ethyl acetate vs φ_ea (%)")
    ax[1, 2].set_xlabel("phi_ea")
    ax[1, 2].set_ylabel("%")
    ax[1, 2].grid(alpha=0.3)
    ax[1, 2].legend(fontsize=8)

plt.suptitle(
    f"Soft coupling ethyl acetate (φ={PHI_ETHYL_ACETATE_STATIC}) + pairwise – snapshot-level",
    fontsize=11, y=1.01,
)
plt.tight_layout()
plt.show()

# ─────────────────────────────────────────────────────────────────────────────────
# G. GUARDADO
# ─────────────────────────────────────────────────────────────────────────────────

if SAVE_EA_CSV:
    _ea_time_file   = OUT_DIR / "ea_soft_coupling_time_comparison_cell28.csv"
    _ea_conc_file   = OUT_DIR / "ea_soft_coupling_concentration_cell28.csv"
    _ea_pool_file   = OUT_DIR / "ea_soft_coupling_pool_distribution_cell28.csv"
    _ea_sens_file   = OUT_DIR / "ea_soft_coupling_sensitivity_cell28.csv"

    ea_comparison_df.to_csv(_ea_time_file, index=False)
    ea_conc_df.to_csv(_ea_conc_file, index=False)
    ea_pool_df.to_csv(_ea_pool_file, index=False)
    if not ea_sensitivity_df.empty:
        ea_sensitivity_df.to_csv(_ea_sens_file, index=False)

    print("\nArchivos guardados:")
    for fp in [_ea_time_file, _ea_conc_file, _ea_pool_file, _ea_sens_file]:
        print(f"  {fp}")

print("\n[Variables disponibles]")
for vn in ["ea_comparison_df", "ea_conc_df", "ea_pool_df", "ea_sensitivity_df", "snapshot_cache_ea"]:
    print(f"  {vn}")

print("\n" + "=" * 100)
print("OK: soft coupling ethyl acetate integrado sobre pairwise (snapshot-level).")
print(f"  PHI_ETHYL_ACETATE_STATIC = {PHI_ETHYL_ACETATE_STATIC}  ← ajusta si el total o la fracción no son razonables.")
print("  Rango recomendado: 0.002 – 0.015")
print("=" * 100)