# ======================================================================================
# MODULO 4 - Funcion objetivo MULTI-EXPERIMENTO (per-observable RMSE-norm + cap)
# ======================================================================================
# Calibra contra todos los ensayos en SELECTED_DATASETS simultaneamente. Para cada
# dataset corre forward_simulate con sus IC + T-profile + pulsos inyectados via
# set_process_inputs (Fases 2-4), interpola la sim en los tiempos experimentales
# y suma residuos ponderados entre todos los ensayos.

OBJ_PENALTY_FAIL = 1e10        # penalty si TODOS los forwards fallan
MIN_SIGMA_FRAC   = 0.05        # sigma minimo relativo si no viene en datos
OBS_J_CAP        = 100.0       # cap por (dataset, observable) en RMSE_norm


def _build_observable_scales(exp_df, sim_nominal):
    """Escala por observable = max abs experimental (sobre TODOS los datasets)."""
    scales = {}
    for v in TARGET_OBSERVABLES:
        sub = exp_df[exp_df["variable"] == v]
        s = float(np.nanmax(np.abs(sub["value"].to_numpy()))) if len(sub) else 1.0
        if (not np.isfinite(s)) or s < 1e-9:
            arr = sim_nominal.get("obs", {}).get(v) if sim_nominal else None
            s = float(np.nanmax(np.abs(arr))) if arr is not None and np.any(np.isfinite(arr)) else 1.0
            if s < 1e-9:
                s = 1.0
        scales[v] = s
    return scales


OBS_SCALES = _build_observable_scales(
    exp_df,
    sim_nom if (isinstance(globals().get("sim_nom"), dict) and sim_nom.get("status") == "ok") else None,
)

print("[OBJ] Escalas por observable (global, multi-ensayo):")
for k, v in OBS_SCALES.items():
    print(f"  {k:24s} scale = {v:.4g}")
print(f"[OBJ] Cap por (dataset, obs) RMSE_norm = {OBS_J_CAP}")


def residuals_single_dataset(sim_result, exp_df_ds):
    """
    Residuos por observable para UN experimento (helper de diagnostico y multi-exp).
    Retorna (dict {var: residuos_norm}, dict {var: RMSE_capped}, J_total_ds).
    """
    if sim_result.get("status") != "ok":
        return None, None, OBS_J_CAP * max(1, len(TARGET_OBSERVABLES))

    t_sim = sim_result["t_h"]
    out_res = {}
    out_rmse = {}

    for v in TARGET_OBSERVABLES:
        sub = exp_df_ds[exp_df_ds["variable"] == v]
        if not len(sub):
            continue
        y_sim_v = sim_result["obs"].get(v)
        if y_sim_v is None or not np.any(np.isfinite(y_sim_v)):
            out_rmse[v] = OBS_J_CAP
            out_res[v]  = np.array([OBS_J_CAP])
            continue
        t_exp = sub["t_h"].to_numpy(dtype=float)
        y_exp = sub["value"].to_numpy(dtype=float)
        sigma = (sub["sigma"].to_numpy(dtype=float)
                 if "sigma" in sub.columns else np.full_like(y_exp, np.nan))
        mask_sim = np.isfinite(y_sim_v)
        if mask_sim.sum() < 2:
            out_rmse[v] = OBS_J_CAP
            out_res[v]  = np.array([OBS_J_CAP])
            continue

        y_interp = np.interp(t_exp, t_sim[mask_sim], y_sim_v[mask_sim])
        sigma_eff = np.where(np.isfinite(sigma) & (sigma > 0), sigma,
                             np.maximum(MIN_SIGMA_FRAC * np.abs(y_exp),
                                        MIN_SIGMA_FRAC * OBS_SCALES[v]))
        scale_v = OBS_SCALES[v]
        res_v = (y_interp - y_exp) / (sigma_eff * scale_v)
        rmse_norm = float(np.sqrt(np.mean(res_v ** 2)))
        rmse_capped = min(rmse_norm, OBS_J_CAP)
        out_res[v]  = res_v
        out_rmse[v] = rmse_capped

    J_ds = 0.0
    for v, r in out_rmse.items():
        J_ds += OBSERVABLE_WEIGHTS.get(v, 1.0) * r
    return out_res, out_rmse, float(J_ds)


# --- Retrocompat: firma (sim_result,) -> (res_dict, J_total) -------------------------
def residuals_by_observable(sim_result, exp_df_local=None):
    """Mantiene API legacy (single-dataset). Usado por post-procesamiento."""
    if exp_df_local is None:
        exp_df_local = exp_df
    res, rmse, J = residuals_single_dataset(sim_result, exp_df_local)
    return res, J


def forward_simulate_multi(theta_unit, datasets=None, *, nfe=None, ncp=None, n_iter=None):
    """
    Corre forward_simulate por cada dataset en `datasets` (default: SELECTED_DATASETS)
    con sus entradas de proceso inyectadas. Retorna dict {ds_id: sim_result}.
    """
    if datasets is None:
        datasets = list(SELECTED_DATASETS)
    sims = {}
    nfe    = NFE_OPT   if nfe    is None else nfe
    ncp    = NCP_OPT   if ncp    is None else ncp
    n_iter = NITER_OPT if n_iter is None else n_iter
    for ds in datasets:
        meta = experiments_meta[ds]
        snap = set_process_inputs(ds)
        try:
            sim = forward_simulate(theta_unit, nfe=nfe, ncp=ncp,
                                   t_end=float(meta["t_end"]),
                                   n_iter=n_iter, silent=True)
        finally:
            restore_process_inputs(snap)
        sim["dataset_id"] = ds
        sims[ds] = sim
    return sims


def objective_multi(theta_unit, *, datasets=None, return_breakdown=False,
                    nfe=None, ncp=None, n_iter=None):
    """
    Objetivo multi-experimento: suma de J por ensayo.
    Si TODOS los forwards fallan -> OBJ_PENALTY_FAIL.
    Si algunos fallan -> esos aportan cap_total (OBS_J_CAP * n_obs_visible).
    """
    if datasets is None:
        datasets = list(SELECTED_DATASETS)
    if not datasets:
        return OBJ_PENALTY_FAIL if not return_breakdown else (OBJ_PENALTY_FAIL, {})

    sims = forward_simulate_multi(theta_unit, datasets=datasets,
                                  nfe=nfe, ncp=ncp, n_iter=n_iter)

    breakdown = {}
    any_ok = False
    J_sum  = 0.0
    for ds in datasets:
        sim = sims[ds]
        sub = exp_df[exp_df["dataset_id"] == ds]
        _, rmse_map, J_ds = residuals_single_dataset(sim, sub)
        if sim.get("status") == "ok":
            any_ok = True
        breakdown[ds] = {
            "status":   sim.get("status"),
            "J":        float(J_ds),
            "rmse":     rmse_map or {},
            "error":    sim.get("error"),
        }
        J_sum += float(J_ds)

    if not any_ok:
        J_sum = OBJ_PENALTY_FAIL

    if return_breakdown:
        return float(J_sum), breakdown
    return float(J_sum)


def objective(theta_unit, *, nfe=NFE_OPT, ncp=NCP_OPT, t_end=TEND_OPT, n_iter=NITER_OPT):
    """Wrapper PSO — multi-experimento (ignora t_end, usa t_end por ensayo)."""
    J = objective_multi(theta_unit, nfe=nfe, ncp=ncp, n_iter=n_iter)
    if not np.isfinite(J):
        return OBJ_PENALTY_FAIL
    return float(J)


# --- Evaluacion en punto nominal -------------------------------------------------------
print(f"\n[OBJ] Evaluando theta_nominal multi-experimento ({len(SELECTED_DATASETS)} ensayos)...")
J_nom, break_nom = objective_multi(theta_nominal_unit, return_breakdown=True)
print(f"[OBJ] J_total(theta_nominal) = {J_nom:.6g}")
print("[OBJ] Descomposicion por ensayo:")
for ds, info in break_nom.items():
    if info["status"] != "ok":
        print(f"  {ds:10s} FAIL ({info['error']})  J_ds={info['J']:.4g}")
        continue
    rmse_str = "  ".join(f"{v[:6]}={r:.3g}" for v, r in info["rmse"].items())
    print(f"  {ds:10s} OK   J_ds={info['J']:8.4g}  | {rmse_str}")
