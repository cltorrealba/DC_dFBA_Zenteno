# Restricción opcional para favorecer acetate esters en snapshots dFBA/FVA
print("\n" + "=" * 100)
print("Acetate esters: restricción fisiológica opcional por snapshots")
print("=" * 100)

required_symbols = [
    "model", "OUT_DIR",
    "states_df_phase", "flux_df_phase", "phase_meta_df",
    "_build_phase_snapshot_model",
    "available_aroma_targets",
    "OBJ_ID", "ATPM_ID", "ETH_ID",
]
missing_symbols = [s for s in required_symbols if s not in globals()]
if missing_symbols:
    raise RuntimeError(f"Faltan símbolos requeridos en sesión: {missing_symbols}")

# --------------------------------------------------------------------------------------
# Configuración
# --------------------------------------------------------------------------------------
ENABLE_ESTER_CONSTRAINT = True
ESTER_CONSTRAINT_MODE = "dynamic_coupling"   # {"dynamic_coupling", "fixed_coupling", "off"}
PHI_ESTER_MAX = 0.12
PHI_ESTER_GRID = [0.00, 0.03, 0.06, 0.12]

MU_LOW_THRESHOLD = 0.03
N_LIMIT_THRESHOLD = 5.0 * float(globals().get("N_TOTAL_DEPLETION_THRESHOLD", 1e-3))
ALCOHOL_ACTIVATION_THRESHOLD = 1e-6
ESTER_BIND_TOL = 1e-7

RUN_EXTENDED_FVA = True
FVA_FRACTIONS_ESTER = (1.0, 0.99)
RUN_PHI_SENSITIVITY = True
SAVE_ESTER_CONSTRAINT_CSV = True

print(f"ENABLE_ESTER_CONSTRAINT = {ENABLE_ESTER_CONSTRAINT}")
print(f"ESTER_CONSTRAINT_MODE   = {ESTER_CONSTRAINT_MODE}")
print(f"PHI_ESTER_MAX          = {PHI_ESTER_MAX}")
print(f"PHI_ESTER_GRID         = {PHI_ESTER_GRID}")

# --------------------------------------------------------------------------------------
# Helpers locales
# --------------------------------------------------------------------------------------
def _slug(txt):
    return re.sub(r"[^0-9a-zA-Z]+", "_", str(txt)).strip("_").lower()

def _rxn_text_local(rxn):
    met_names = " ".join(getattr(m, "name", "") for m in rxn.metabolites)
    return f"{rxn.id} {rxn.name or ''} {met_names}".lower()

def _is_exchange_like_local(rxn):
    return ("exchange" in (rxn.name or "").lower()) or (len(rxn.metabolites) == 1)

def _find_best_rxn_id(keyword_sets, manual_id=None, exchange_preferred=True):
    rxn_ids_all = {r.id for r in model.reactions}
    if manual_id is not None and manual_id in rxn_ids_all:
        rxn = model.reactions.get_by_id(manual_id)
        return manual_id, [{"score": 0, "rid": rxn.id, "name": rxn.name, "matched": "manual"}]

    hits = []
    for rxn in model.reactions:
        txt = _rxn_text_local(rxn)
        for kws in keyword_sets:
            kws_l = [str(k).lower() for k in kws]
            if all(k in txt for k in kws_l):
                score = len(rxn.metabolites)
                if exchange_preferred and _is_exchange_like_local(rxn):
                    score -= 10
                hits.append(
                    {
                        "score": score,
                        "rid": rxn.id,
                        "name": rxn.name,
                        "matched": " & ".join(kws_l),
                    }
                )
                break
    if not hits:
        return None, []
    hits = sorted(hits, key=lambda x: (x["score"], x["rid"]))
    return hits[0]["rid"], hits[:10]

def _get_exchange_mw_g_per_mmol_local(mdl, rid):
    rxn = mdl.reactions.get_by_id(rid)
    mets = list(rxn.metabolites.keys())
    if len(mets) != 1:
        return np.nan
    met = mets[0]
    try:
        return float(met.formula_weight) / 1000.0
    except Exception:
        return np.nan

def _sum_fluxes(sol, rid_list, positive_only=False):
    vals = []
    for rid in rid_list:
        v = float(sol.fluxes.get(rid, 0.0))
        vals.append(max(0.0, v) if positive_only else v)
    return float(np.sum(vals))

def _safe_flux(sol, rid):
    return float(sol.fluxes.get(rid, np.nan)) if getattr(sol, "status", "") == "optimal" else np.nan

# --------------------------------------------------------------------------------------
# A) Inspección y mapeo de IDs relevantes
# --------------------------------------------------------------------------------------
rxn_ids_all = {r.id for r in model.reactions}

# 1) Acetate esters por exchanges ya detectados aguas arriba
acetate_ester_targets = {
    rid: label
    for rid, label in available_aroma_targets.items()
    if "acetate" in str(label).lower()
}

if not acetate_ester_targets:
    raise RuntimeError("No se detectaron acetate esters en available_aroma_targets.")

# 2) Búsqueda de alcoholes superiores relacionados
found_aromas_local = globals().get("found_aromas", {})
manual_alcohol_map = {
    "ethyl acetate exchange": ETH_ID if ETH_ID in rxn_ids_all else None,
    "isobutyl acetate exchange": found_aromas_local.get("isobutanol", None),
    "isoamyl acetate exchange": found_aromas_local.get("isoamyl", None),
    "phenethyl acetate exchange": found_aromas_local.get("pea", None),
    "hexyl acetate exchange": None,
    "benzyl-acetate exchange": None,
}

alcohol_search_map = {
    "ethyl acetate exchange": [["ethanol", "exchange"]],
    "isobutyl acetate exchange": [["isobutanol", "exchange"], ["2-methyl-1-propanol", "exchange"]],
    "isoamyl acetate exchange": [["isoamyl alcohol", "exchange"], ["3-methyl-1-butanol", "exchange"], ["isoamyl", "exchange"]],
    "phenethyl acetate exchange": [["2-phenylethanol", "exchange"], ["phenylethanol", "exchange"], ["pea", "exchange"]],
    "hexyl acetate exchange": [["hexanol", "exchange"], ["hexan-1-ol", "exchange"]],
    "benzyl-acetate exchange": [["benzyl alcohol", "exchange"], ["benzylalcohol", "exchange"]],
}

ester_alcohol_rows = []
ester_to_alcohol = {}
for ester_rid, ester_label in acetate_ester_targets.items():
    alcohol_rid, alcohol_hits = _find_best_rxn_id(
        alcohol_search_map.get(ester_label, []),
        manual_id=manual_alcohol_map.get(ester_label),
        exchange_preferred=True,
    )
    ester_to_alcohol[ester_rid] = alcohol_rid
    ester_alcohol_rows.append(
        {
            "ester_rid": ester_rid,
            "ester_label": ester_label,
            "alcohol_rid": alcohol_rid,
            "alcohol_name": model.reactions.get_by_id(alcohol_rid).name if alcohol_rid in rxn_ids_all else None,
            "n_hits": len(alcohol_hits),
            "top_match": alcohol_hits[0]["matched"] if alcohol_hits else None,
        }
    )

ester_alcohol_map_df = pd.DataFrame(ester_alcohol_rows)
ester_ids = list(acetate_ester_targets.keys())
higher_alcohol_ids = list(dict.fromkeys([rid for rid in ester_to_alcohol.values() if rid is not None and rid in rxn_ids_all]))

# 3) Metabolitos proxy: acetyl-CoA / CoA
proxy_met_rows = []
for met in model.metabolites:
    nm = (met.name or "").lower()
    if ("acetyl-coa" in nm) or (nm == "coa") or ("coenzyme a" in nm) or ("acetyl coenzyme a" in nm):
        proxy_met_rows.append(
            {
                "met_id": met.id,
                "name": met.name,
                "formula": getattr(met, "formula", None),
                "n_reactions": len(getattr(met, "reactions", [])),
            }
        )
proxy_metabolites_df = pd.DataFrame(proxy_met_rows).sort_values(["name", "met_id"], kind="stable") if proxy_met_rows else pd.DataFrame()

# 4) Reacciones internas candidatas ligadas a esterificación/acetyl-CoA
proxy_rxn_rows = []
proxy_rxn_keywords = ["acetyltransferase", "ester", "acetyl-coa", "acetyl coa", "alcohol acetyltransferase", "coa transferase"]
for rxn in model.reactions:
    txt = _rxn_text_local(rxn)
    matched = [kw for kw in proxy_rxn_keywords if kw in txt]
    if matched:
        proxy_rxn_rows.append(
            {
                "rid": rxn.id,
                "name": rxn.name,
                "n_mets": len(rxn.metabolites),
                "matched_keywords": ", ".join(matched),
            }
        )
proxy_reactions_df = pd.DataFrame(proxy_rxn_rows).sort_values(["rid"], kind="stable") if proxy_rxn_rows else pd.DataFrame()

print("\n[Mapeo ester ↔ alcohol relacionado]")
display(ester_alcohol_map_df)

if not proxy_metabolites_df.empty:
    print("\n[Metabolitos proxy detectados: acetyl-CoA / CoA]")
    display(proxy_metabolites_df)

if not proxy_reactions_df.empty:
    print("\n[Reacciones candidatas ligadas a esterificación / acetyl-CoA]")
    display(proxy_reactions_df.head(20))

if not higher_alcohol_ids:
    print("[WARN] No se detectaron alcoholes relacionados robustos; la restricción quedará en fallback inactivo.")

# --------------------------------------------------------------------------------------
# B) Wrapper incremental sobre _build_phase_snapshot_model(...)
# --------------------------------------------------------------------------------------
def _compute_phi_ester(state_row, sol_ref, mode, alcohol_sum_pos):
    if (not ENABLE_ESTER_CONSTRAINT) or (ESTER_CONSTRAINT_MODE == "off"):
        return 0.0, 0.0, {
            "turnover": 0.0,
            "n_limited": 0.0,
            "low_growth": 0.0,
            "anaerobic": 0.0,
            "alcohol_present": 0.0,
        }

    n_val = float(state_row["N_gN_L"]) if "N_gN_L" in state_row.index else np.nan
    o2_val = float(state_row["O2_g_L"]) if "O2_g_L" in state_row.index else 0.0
    mu_val = float(sol_ref.fluxes.get(OBJ_ID, 0.0)) if getattr(sol_ref, "status", "") == "optimal" else np.nan

    flags = {
        "turnover": float(mode == "TURNOVER_ATPM"),
        "n_limited": float(np.isfinite(n_val) and n_val <= N_LIMIT_THRESHOLD),
        "low_growth": float(np.isfinite(mu_val) and mu_val <= MU_LOW_THRESHOLD),
        "anaerobic": float(np.isfinite(o2_val) and o2_val <= max(1e-9, 10.0 * float(globals().get("O2_DEPLETION_THRESHOLD", 1e-6)))),
        "alcohol_present": float(alcohol_sum_pos >= ALCOHOL_ACTIVATION_THRESHOLD),
    }

    proxy_score = (
        0.40 * flags["turnover"]
        + 0.25 * flags["n_limited"]
        + 0.20 * flags["low_growth"]
        + 0.10 * flags["anaerobic"]
        + 0.05 * flags["alcohol_present"]
    )
    proxy_score = float(min(1.0, max(0.0, proxy_score)))

    if ESTER_CONSTRAINT_MODE == "fixed_coupling":
        phi = PHI_ESTER_MAX if flags["alcohol_present"] else 0.0
    else:
        phi = PHI_ESTER_MAX * proxy_score
        if not flags["alcohol_present"]:
            phi *= 0.25

    return float(phi), float(proxy_score), flags

def _add_ester_coupling_constraint(mtmp, phi_ester):
    if phi_ester <= 0.0 or not ester_ids or not higher_alcohol_ids:
        return None, False

    rb = {r.id: r for r in mtmp.reactions}
    ester_expr = sum(rb[rid].flux_expression for rid in ester_ids if rid in rb)
    alcohol_expr = sum(rb[rid].flux_expression for rid in higher_alcohol_ids if rid in rb)

    cname = f"ester_coupling_phi_{str(phi_ester).replace('.', '_')}"
    cons = mtmp.problem.Constraint(ester_expr - float(phi_ester) * alcohol_expr, lb=0.0, name=cname)
    mtmp.add_cons_vars([cons])
    return cons, True

def _build_phase_snapshot_model_with_ester_constraint(t_phase, *, enable_constraint=True, phi_override=None):
    snap = _build_phase_snapshot_model(float(t_phase))
    m_base = snap["model"]
    sol_base = snap["solution"]
    state_row = snap["state_row"]
    mode = snap["mode"]
    objective_id = snap["objective_id"]

    alcohol_sum_base_pos = _sum_fluxes(sol_base, higher_alcohol_ids, positive_only=True) if getattr(sol_base, "status", "") == "optimal" else np.nan
    ester_sum_base_pos = _sum_fluxes(sol_base, ester_ids, positive_only=True) if getattr(sol_base, "status", "") == "optimal" else np.nan

    phi_auto, proxy_score, proxy_flags = _compute_phi_ester(state_row, sol_base, mode, alcohol_sum_base_pos if np.isfinite(alcohol_sum_base_pos) else 0.0)
    phi_ester = float(phi_auto if phi_override is None else phi_override)

    constraint_active = bool(
        enable_constraint
        and ESTER_CONSTRAINT_MODE != "off"
        and phi_ester > 0.0
        and len(ester_ids) > 0
        and len(higher_alcohol_ids) > 0
    )

    m_con = m_base.copy()
    sol_con = sol_base
    lhs_margin = np.nan
    constraint_binding = False
    solved_constrained = False

    if constraint_active:
        _cons, ok = _add_ester_coupling_constraint(m_con, phi_ester)
        if ok:
            try:
                sol_con = m_con.optimize()
                solved_constrained = True
            except Exception:
                pass

    if getattr(sol_con, "status", "") == "optimal":
        ester_sum_con_pos = _sum_fluxes(sol_con, ester_ids, positive_only=True)
        alcohol_sum_con_pos = _sum_fluxes(sol_con, higher_alcohol_ids, positive_only=True)
        lhs_margin = _sum_fluxes(sol_con, ester_ids, positive_only=False) - phi_ester * _sum_fluxes(sol_con, higher_alcohol_ids, positive_only=False)
        constraint_binding = bool(constraint_active and np.isfinite(lhs_margin) and lhs_margin <= ESTER_BIND_TOL)
    else:
        ester_sum_con_pos = np.nan
        alcohol_sum_con_pos = np.nan

    return {
        "t_h": float(t_phase),
        "mode": mode,
        "objective_id": objective_id,
        "state_row": state_row,
        "model_baseline": m_base,
        "model_constrained": m_con,
        "solution_baseline": sol_base,
        "solution_constrained": sol_con,
        "phi_ester": phi_ester,
        "proxy_score": proxy_score,
        "proxy_flags": proxy_flags,
        "constraint_active": constraint_active,
        "constraint_binding": constraint_binding,
        "lhs_margin": lhs_margin,
        "ester_sum_base_pos": ester_sum_base_pos,
        "alcohol_sum_base_pos": alcohol_sum_base_pos,
        "ester_sum_con_pos": ester_sum_con_pos,
        "alcohol_sum_con_pos": alcohol_sum_con_pos,
        "solved_constrained": solved_constrained,
    }

# --------------------------------------------------------------------------------------
# C) Comparación baseline vs constrained en el tiempo
# --------------------------------------------------------------------------------------
analysis_times = flux_df_phase["t_h"].to_numpy(dtype=float)
state_times = states_df_phase["t_h"].to_numpy(dtype=float)
phase_name_by_time = {float(r["t_h"]): str(r["fase"]) for _, r in phase_meta_df.iterrows()}

ester_mw_map = {rid: _get_exchange_mw_g_per_mmol_local(model, rid) for rid in ester_ids}
higher_alcohol_label_map = {rid: model.reactions.get_by_id(rid).name for rid in higher_alcohol_ids}

ester_conc_baseline = {rid: np.zeros(len(state_times), dtype=float) for rid in ester_ids}
ester_conc_constrained = {rid: np.zeros(len(state_times), dtype=float) for rid in ester_ids}

comparison_rows = []
snapshot_cache_baseline = {}
snapshot_cache_constrained = {}

for k, t_k in enumerate(analysis_times):
    res_cmp = _build_phase_snapshot_model_with_ester_constraint(
        t_k,
        enable_constraint=ENABLE_ESTER_CONSTRAINT,
        phi_override=None,
    )

    snapshot_cache_baseline[float(t_k)] = {
        "model": res_cmp["model_baseline"],
        "solution": res_cmp["solution_baseline"],
    }
    snapshot_cache_constrained[float(t_k)] = {
        "model": res_cmp["model_constrained"],
        "solution": res_cmp["solution_constrained"],
    }

    row_state = res_cmp["state_row"]
    sol_b = res_cmp["solution_baseline"]
    sol_c = res_cmp["solution_constrained"]

    obj_id = res_cmp["objective_id"]
    mu_b = _safe_flux(sol_b, OBJ_ID)
    mu_c = _safe_flux(sol_c, OBJ_ID)
    atpm_b = _safe_flux(sol_b, ATPM_ID)
    atpm_c = _safe_flux(sol_c, ATPM_ID)
    obj_b = _safe_flux(sol_b, obj_id)
    obj_c = _safe_flux(sol_c, obj_id)

    out_row = {
        "t_h": float(t_k),
        "fase": phase_name_by_time.get(float(t_k), None),
        "mode": res_cmp["mode"],
        "objective_id": obj_id,
        "objective_value_baseline": obj_b,
        "objective_value_constrained": obj_c,
        "mu_h_inv_baseline": mu_b,
        "mu_h_inv_constrained": mu_c,
        "v_ATPM_baseline": atpm_b,
        "v_ATPM_constrained": atpm_c,
        "ester_total_baseline": res_cmp["ester_sum_base_pos"],
        "ester_total_constrained": res_cmp["ester_sum_con_pos"],
        "higher_alcohol_total_baseline": res_cmp["alcohol_sum_base_pos"],
        "higher_alcohol_total_constrained": res_cmp["alcohol_sum_con_pos"],
        "proxy_score": res_cmp["proxy_score"],
        "phi_ester": res_cmp["phi_ester"],
        "constraint_active": res_cmp["constraint_active"],
        "constraint_binding": res_cmp["constraint_binding"],
        "lhs_margin": res_cmp["lhs_margin"],
        "delta_objective_vs_baseline": obj_c - obj_b if np.isfinite(obj_b) and np.isfinite(obj_c) else np.nan,
        "delta_mu_vs_baseline": mu_c - mu_b if np.isfinite(mu_b) and np.isfinite(mu_c) else np.nan,
        "delta_ester_total_vs_baseline": res_cmp["ester_sum_con_pos"] - res_cmp["ester_sum_base_pos"] if np.isfinite(res_cmp["ester_sum_con_pos"]) and np.isfinite(res_cmp["ester_sum_base_pos"]) else np.nan,
        "proxy_turnover": res_cmp["proxy_flags"]["turnover"],
        "proxy_n_limited": res_cmp["proxy_flags"]["n_limited"],
        "proxy_low_growth": res_cmp["proxy_flags"]["low_growth"],
        "proxy_anaerobic": res_cmp["proxy_flags"]["anaerobic"],
        "proxy_alcohol_present": res_cmp["proxy_flags"]["alcohol_present"],
        "N_gN_L": float(row_state["N_gN_L"]) if "N_gN_L" in row_state.index else np.nan,
        "E_g_L": float(row_state["E_g_L"]) if "E_g_L" in row_state.index else np.nan,
        "O2_g_L": float(row_state["O2_g_L"]) if "O2_g_L" in row_state.index else np.nan,
        "X_gDW_L": float(row_state["X_gDW_L"]) if "X_gDW_L" in row_state.index else np.nan,
    }

    for rid, label in acetate_ester_targets.items():
        out_row[f"ester_{rid}_baseline"] = _safe_flux(sol_b, rid)
        out_row[f"ester_{rid}_constrained"] = _safe_flux(sol_c, rid)
        out_row[f"ester_{rid}_delta"] = out_row[f"ester_{rid}_constrained"] - out_row[f"ester_{rid}_baseline"] if np.isfinite(out_row[f"ester_{rid}_baseline"]) and np.isfinite(out_row[f"ester_{rid}_constrained"]) else np.nan

    for rid in higher_alcohol_ids:
        out_row[f"alcohol_{rid}_baseline"] = _safe_flux(sol_b, rid)
        out_row[f"alcohol_{rid}_constrained"] = _safe_flux(sol_c, rid)
        out_row[f"alcohol_{rid}_delta"] = out_row[f"alcohol_{rid}_constrained"] - out_row[f"alcohol_{rid}_baseline"] if np.isfinite(out_row[f"alcohol_{rid}_baseline"]) and np.isfinite(out_row[f"alcohol_{rid}_constrained"]) else np.nan

    comparison_rows.append(out_row)

    # integración macroscópica acumulada de acetate esters
    if k + 1 < len(state_times):
        dt_k = float(state_times[k + 1] - state_times[k])
        Xk = float(row_state["X_gDW_L"]) if "X_gDW_L" in row_state.index else 0.0

        for rid in ester_ids:
            mw = ester_mw_map.get(rid, np.nan)
            vb = max(0.0, _safe_flux(sol_b, rid)) if np.isfinite(_safe_flux(sol_b, rid)) else 0.0
            vc = max(0.0, _safe_flux(sol_c, rid)) if np.isfinite(_safe_flux(sol_c, rid)) else 0.0

            ester_conc_baseline[rid][k + 1] = ester_conc_baseline[rid][k]
            ester_conc_constrained[rid][k + 1] = ester_conc_constrained[rid][k]

            if np.isfinite(mw):
                ester_conc_baseline[rid][k + 1] += 1000.0 * mw * vb * Xk * dt_k
                ester_conc_constrained[rid][k + 1] += 1000.0 * mw * vc * Xk * dt_k

ester_constraint_time_df = pd.DataFrame(comparison_rows)

ester_constraint_conc_df = pd.DataFrame({"t_h": state_times})
for rid, label in acetate_ester_targets.items():
    sfx = _slug(label)
    ester_constraint_conc_df[f"{sfx}_baseline_mg_L"] = ester_conc_baseline[rid]
    ester_constraint_conc_df[f"{sfx}_constrained_mg_L"] = ester_conc_constrained[rid]
    ester_constraint_conc_df[f"{sfx}_delta_mg_L"] = ester_conc_constrained[rid] - ester_conc_baseline[rid]

baseline_conc_cols = [c for c in ester_constraint_conc_df.columns if c.endswith("_baseline_mg_L")]
constrained_conc_cols = [c for c in ester_constraint_conc_df.columns if c.endswith("_constrained_mg_L")]
ester_constraint_conc_df["total_acetate_esters_baseline_mg_L"] = ester_constraint_conc_df[baseline_conc_cols].sum(axis=1)
ester_constraint_conc_df["total_acetate_esters_constrained_mg_L"] = ester_constraint_conc_df[constrained_conc_cols].sum(axis=1)
ester_constraint_conc_df["total_acetate_esters_delta_mg_L"] = (
    ester_constraint_conc_df["total_acetate_esters_constrained_mg_L"]
    - ester_constraint_conc_df["total_acetate_esters_baseline_mg_L"]
)

phase_times = set(float(x) for x in phase_meta_df["t_h"].tolist())
ester_constraint_phase_df = ester_constraint_time_df[ester_constraint_time_df["t_h"].isin(phase_times)].copy()

print("\n[Comparativa por fase: baseline vs constrained]")
display(ester_constraint_phase_df)

# --------------------------------------------------------------------------------------
# D1) Gráficos nuevos
# --------------------------------------------------------------------------------------
fig, ax = plt.subplots(2, 2, figsize=(14, 8), sharex=True)

ax[0, 0].plot(
    ester_constraint_time_df["t_h"],
    ester_constraint_time_df["ester_total_baseline"],
    lw=2,
    label="baseline",
)
ax[0, 0].plot(
    ester_constraint_time_df["t_h"],
    ester_constraint_time_df["ester_total_constrained"],
    lw=2,
    ls="--",
    label="constrained",
)
ax[0, 0].set_title("Flujo total acetate esters")
ax[0, 0].set_ylabel("mmol/gDW/h")
ax[0, 0].grid(alpha=0.3)
ax[0, 0].legend()

ax[0, 1].plot(
    ester_constraint_conc_df["t_h"],
    ester_constraint_conc_df["total_acetate_esters_baseline_mg_L"],
    lw=2,
    label="baseline",
)
ax[0, 1].plot(
    ester_constraint_conc_df["t_h"],
    ester_constraint_conc_df["total_acetate_esters_constrained_mg_L"],
    lw=2,
    ls="--",
    label="constrained",
)
ax[0, 1].set_title("Acumulado total acetate esters")
ax[0, 1].set_ylabel("mg/L")
ax[0, 1].grid(alpha=0.3)
ax[0, 1].legend()

ax[1, 0].plot(
    ester_constraint_time_df["t_h"],
    ester_constraint_time_df["proxy_score"],
    lw=2,
    label="proxy_score",
)
ax[1, 0].plot(
    ester_constraint_time_df["t_h"],
    ester_constraint_time_df["phi_ester"],
    lw=2,
    ls="--",
    label="phi_ester",
)
ax[1, 0].set_title("Proxy fisiológica y phi_ester")
ax[1, 0].set_xlabel("t (h)")
ax[1, 0].grid(alpha=0.3)
ax[1, 0].legend()

ax[1, 1].plot(
    ester_constraint_time_df["t_h"],
    ester_constraint_time_df["mu_h_inv_baseline"],
    lw=2,
    label="mu baseline",
)
ax[1, 1].plot(
    ester_constraint_time_df["t_h"],
    ester_constraint_time_df["mu_h_inv_constrained"],
    lw=2,
    ls="--",
    label="mu constrained",
)
ax2 = ax[1, 1].twinx()
ax2.plot(
    ester_constraint_time_df["t_h"],
    ester_constraint_time_df["v_ATPM_baseline"],
    lw=1.6,
    color="tab:orange",
    alpha=0.8,
    label="ATPM baseline",
)
ax2.plot(
    ester_constraint_time_df["t_h"],
    ester_constraint_time_df["v_ATPM_constrained"],
    lw=1.6,
    ls="--",
    color="tab:red",
    alpha=0.8,
    label="ATPM constrained",
)
ax[1, 1].set_title("Trade-off con crecimiento / ATPM")
ax[1, 1].set_xlabel("t (h)")
ax[1, 1].set_ylabel("mu (h^-1)")
ax2.set_ylabel("ATPM (mmol/gDW/h)")
ax[1, 1].grid(alpha=0.3)

lines1, labels1 = ax[1, 1].get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax[1, 1].legend(lines1 + lines2, labels1 + labels2, loc="best")

plt.tight_layout()
plt.show()

# Figura adicional: principales acetate esters individuales
peak_by_ester = {}
for rid in ester_ids:
    col = f"ester_{rid}_constrained"
    peak_by_ester[rid] = float(np.nanmax(ester_constraint_time_df[col].to_numpy(dtype=float))) if col in ester_constraint_time_df else 0.0

top_ester_ids = sorted(peak_by_ester, key=peak_by_ester.get, reverse=True)[: min(4, len(ester_ids))]

if top_ester_ids:
    ncols = min(2, len(top_ester_ids))
    nrows = int(np.ceil(len(top_ester_ids) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.2 * ncols, 3.8 * nrows), sharex=True)
    axes = np.atleast_1d(axes).ravel()

    for i, rid in enumerate(top_ester_ids):
        label = acetate_ester_targets[rid]
        axes[i].plot(ester_constraint_time_df["t_h"], ester_constraint_time_df[f"ester_{rid}_baseline"], lw=2, label="baseline")
        axes[i].plot(ester_constraint_time_df["t_h"], ester_constraint_time_df[f"ester_{rid}_constrained"], lw=2, ls="--", label="constrained")
        axes[i].set_title(f"{label} ({rid})")
        axes[i].set_xlabel("t (h)")
        axes[i].set_ylabel("mmol/gDW/h")
        axes[i].grid(alpha=0.3)
        axes[i].legend()

    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    plt.tight_layout()
    plt.show()

# --------------------------------------------------------------------------------------
# D2) FVA ampliado: baseline vs constrained
# --------------------------------------------------------------------------------------
ester_constraint_fva_df = pd.DataFrame()

if RUN_EXTENDED_FVA:
    fva_rows = []

    for _, phase_row in phase_meta_df.iterrows():
        t_phase = float(phase_row["t_h"])
        fase = str(phase_row["fase"])

        snap_b = _build_phase_snapshot_model(t_phase)
        snap_c = _build_phase_snapshot_model_with_ester_constraint(
            t_phase,
            enable_constraint=ENABLE_ESTER_CONSTRAINT,
            phi_override=None,
        )

        sol_b = snap_b["solution"]
        sol_c = snap_c["solution_constrained"]

        for frac in FVA_FRACTIONS_ESTER:
            # baseline
            if getattr(sol_b, "status", "") == "optimal":
                try:
                    fva_b = cobra.flux_analysis.flux_variability_analysis(
                        snap_b["model"],
                        reaction_list=ester_ids,
                        fraction_of_optimum=float(frac),
                    )
                    for rid in ester_ids:
                        fva_rows.append(
                            {
                                "scenario": "baseline",
                                "fase": fase,
                                "t_h": t_phase,
                                "fraction_of_optimum": float(frac),
                                "rid": rid,
                                "label": acetate_ester_targets[rid],
                                "flux_ref": float(sol_b.fluxes.get(rid, np.nan)),
                                "minimum": float(fva_b.loc[rid, "minimum"]),
                                "maximum": float(fva_b.loc[rid, "maximum"]),
                                "range": float(fva_b.loc[rid, "maximum"] - fva_b.loc[rid, "minimum"]),
                            }
                        )
                except Exception:
                    pass

            # constrained
            if getattr(sol_c, "status", "") == "optimal":
                try:
                    fva_c = cobra.flux_analysis.flux_variability_analysis(
                        snap_c["model_constrained"],
                        reaction_list=ester_ids,
                        fraction_of_optimum=float(frac),
                    )
                    for rid in ester_ids:
                        fva_rows.append(
                            {
                                "scenario": "constrained",
                                "fase": fase,
                                "t_h": t_phase,
                                "fraction_of_optimum": float(frac),
                                "rid": rid,
                                "label": acetate_ester_targets[rid],
                                "flux_ref": float(sol_c.fluxes.get(rid, np.nan)),
                                "minimum": float(fva_c.loc[rid, "minimum"]),
                                "maximum": float(fva_c.loc[rid, "maximum"]),
                                "range": float(fva_c.loc[rid, "maximum"] - fva_c.loc[rid, "minimum"]),
                            }
                        )
                except Exception:
                    pass

    ester_constraint_fva_df = pd.DataFrame(fva_rows)

    if not ester_constraint_fva_df.empty:
        print("\n[FVA ampliado: baseline vs constrained]")
        display(ester_constraint_fva_df)

        fva_compare = ester_constraint_fva_df.pivot_table(
            index=["fase", "fraction_of_optimum", "rid", "label"],
            columns="scenario",
            values=["flux_ref", "minimum", "maximum", "range"],
            aggfunc="first",
        )
        print("\n[Resumen FVA pivot]")
        display(fva_compare)

# --------------------------------------------------------------------------------------
# D3) Sensibilidad mínima en phi_ester
# --------------------------------------------------------------------------------------
ester_constraint_sensitivity_df = pd.DataFrame()

if RUN_PHI_SENSITIVITY:
    sens_rows = []
    phi_grid = sorted(set(float(x) for x in PHI_ESTER_GRID if float(x) >= 0.0))

    for _, phase_row in phase_meta_df.iterrows():
        t_phase = float(phase_row["t_h"])
        fase = str(phase_row["fase"])

        for phi_val in phi_grid:
            snap_phi = _build_phase_snapshot_model_with_ester_constraint(
                t_phase,
                enable_constraint=(phi_val > 0.0),
                phi_override=phi_val,
            )
            sol_phi = snap_phi["solution_constrained"]
            status_phi = getattr(sol_phi, "status", "unknown")

            sens_rows.append(
                {
                    "fase": fase,
                    "t_h": t_phase,
                    "phi_ester": phi_val,
                    "status": status_phi,
                    "constraint_active": snap_phi["constraint_active"],
                    "constraint_binding": snap_phi["constraint_binding"],
                    "objective_id": snap_phi["objective_id"],
                    "objective_value": _safe_flux(sol_phi, snap_phi["objective_id"]),
                    "mu_h_inv": _safe_flux(sol_phi, OBJ_ID),
                    "v_ATPM": _safe_flux(sol_phi, ATPM_ID),
                    "ester_total": _sum_fluxes(sol_phi, ester_ids, positive_only=True) if status_phi == "optimal" else np.nan,
                    "higher_alcohol_total": _sum_fluxes(sol_phi, higher_alcohol_ids, positive_only=True) if status_phi == "optimal" else np.nan,
                    "lhs_margin": snap_phi["lhs_margin"],
                }
            )

    ester_constraint_sensitivity_df = pd.DataFrame(sens_rows)
    print("\n[Sensibilidad mínima en phi_ester]")
    display(ester_constraint_sensitivity_df)

# --------------------------------------------------------------------------------------
# D4) Resumen final y guardado
# --------------------------------------------------------------------------------------
summary_rows_ester = []
for rid, label in acetate_ester_targets.items():
    sfx = _slug(label)
    final_base = float(ester_constraint_conc_df[f"{sfx}_baseline_mg_L"].iloc[-1])
    final_con = float(ester_constraint_conc_df[f"{sfx}_constrained_mg_L"].iloc[-1])
    summary_rows_ester.append(
        {
            "rid": rid,
            "label": label,
            "final_baseline_mg_L": final_base,
            "final_constrained_mg_L": final_con,
            "delta_mg_L": final_con - final_base,
        }
    )

ester_constraint_summary_df = pd.DataFrame(summary_rows_ester).sort_values("delta_mg_L", ascending=False, kind="stable")

print("\n[Resumen final acetate esters]")
display(ester_constraint_summary_df)

if SAVE_ESTER_CONSTRAINT_CSV:
    ester_map_file = OUT_DIR / "ester_constraint_mapping_cell24.csv"
    ester_time_file = OUT_DIR / "ester_constraint_time_comparison_cell24.csv"
    ester_conc_file = OUT_DIR / "ester_constraint_concentration_comparison_cell24.csv"
    ester_summary_file = OUT_DIR / "ester_constraint_summary_cell24.csv"

    ester_alcohol_map_df.to_csv(ester_map_file, index=False)
    ester_constraint_time_df.to_csv(ester_time_file, index=False)
    ester_constraint_conc_df.to_csv(ester_conc_file, index=False)
    ester_constraint_summary_df.to_csv(ester_summary_file, index=False)

    print("\nArchivos guardados:")
    print(f"- {ester_map_file}")
    print(f"- {ester_time_file}")
    print(f"- {ester_conc_file}")
    print(f"- {ester_summary_file}")

    if not ester_constraint_fva_df.empty:
        ester_fva_file = OUT_DIR / "ester_constraint_fva_cell24.csv"
        ester_constraint_fva_df.to_csv(ester_fva_file, index=False)
        print(f"- {ester_fva_file}")

    if not ester_constraint_sensitivity_df.empty:
        ester_sens_file = OUT_DIR / "ester_constraint_sensitivity_cell24.csv"
        ester_constraint_sensitivity_df.to_csv(ester_sens_file, index=False)
        print(f"- {ester_sens_file}")

print("\nVariables nuevas:")
print("- ester_alcohol_map_df")
print("- proxy_metabolites_df")
print("- proxy_reactions_df")
print("- ester_constraint_time_df")
print("- ester_constraint_phase_df")
print("- ester_constraint_conc_df")
print("- ester_constraint_summary_df")
print("- ester_constraint_fva_df")
print("- ester_constraint_sensitivity_df")

if ENABLE_ESTER_CONSTRAINT and ESTER_CONSTRAINT_MODE != "off":
    print("\nNota:")
    print("- La restricción aplicada es un acoplamiento lineal sobre exchanges.")
    print("- No fija un flujo absoluto; fuerza consistencia relativa ester/alcohol.")
    print("- Si no hay alcoholes mapeados o phi_ester≈0, el comportamiento queda equivalente al baseline.")