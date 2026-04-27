# ======================================================================================
# MODULO 3d - Verificacion pre-PSO: simulacion nominal por experimento vs datos reales
# ======================================================================================
# Objetivo: antes de lanzar PSO, mostrar overlay de forward_simulate(theta_nominal)
# corrido por-experimento (usando IC reales via set_ic_from_experiment) contra los
# datos de exp_df.
#
# Estado actual (Fases 2-4 implementadas):
#   [x] IC por experimento ruteadas al solver (globales X0_init, G0_init, ...).
#   [x] T-profile por punto de colocacion (via T_PROFILE_DF).
#   [x] Pulsos de N al final del FE donde caen (via PULSE_SCHEDULE).
#
# Si exp_df esta vacio o experiments_meta no tiene ensayos, la celda no hace nada.

import matplotlib.pyplot as plt


def forward_simulate_experiment(theta_unit, dataset_id, *, t_end=None, **fwd_kw):
    """
    Corre forward_simulate con las IC del ensayo `dataset_id` inyectadas como globales.
    t_end por defecto = experiments_meta[ds]['t_end'].
    """
    if dataset_id not in experiments_meta:
        raise KeyError(f"dataset {dataset_id!r} no esta en experiments_meta")
    meta = experiments_meta[dataset_id]
    if t_end is None:
        t_end = float(meta["t_end"])
    snap = set_process_inputs(dataset_id)
    try:
        sim = forward_simulate(theta_unit, t_end=t_end, silent=True, **fwd_kw)
    finally:
        restore_process_inputs(snap)
    sim["dataset_id"] = dataset_id
    return sim


def _residuals_for_dataset(sim, ds_id):
    """RMSE_norm por observable para un experimento (usa OBS_SCALES si existe, sino 1)."""
    if sim.get("status") != "ok":
        return None
    sub_all = exp_df[exp_df["dataset_id"] == ds_id]
    t_sim = sim["t_h"]
    out = {}
    for v in TARGET_OBSERVABLES:
        sub = sub_all[sub_all["variable"] == v]
        if not len(sub):
            continue
        y_sim = sim["obs"].get(v)
        if y_sim is None or not np.any(np.isfinite(y_sim)):
            continue
        t_exp = sub["t_h"].to_numpy(dtype=float)
        y_exp = sub["value"].to_numpy(dtype=float)
        sigma = sub["sigma"].to_numpy(dtype=float)
        mask  = np.isfinite(y_sim)
        if mask.sum() < 2:
            continue
        y_int = np.interp(t_exp, t_sim[mask], y_sim[mask])
        scale = max(abs(y_exp).max(), 1e-9)
        res   = (y_int - y_exp) / (np.where(sigma > 0, sigma, 1e-6) * scale)
        out[v] = float(np.sqrt(np.mean(res ** 2)))
    return out


def plot_nominal_vs_experiment(dataset_id, sim=None, theta_unit=None):
    """Grafico de overlay sim nominal vs datos del ensayo `dataset_id`."""
    if sim is None:
        if theta_unit is None:
            theta_unit = theta_nominal_unit
        sim = forward_simulate_experiment(theta_unit, dataset_id)

    sub_all = exp_df[exp_df["dataset_id"] == dataset_id]
    vars_plot = [v for v in TARGET_OBSERVABLES if v in sub_all["variable"].unique()]
    if not vars_plot:
        print(f"[VERIF] {dataset_id}: sin observables en exp_df, skip.")
        return sim

    ncols = min(4, len(vars_plot))
    nrows = int(np.ceil(len(vars_plot) / ncols))
    fig, axs = plt.subplots(nrows, ncols, figsize=(3.6 * ncols, 2.8 * nrows), dpi=110)
    axs = np.atleast_2d(axs).ravel()

    ok = sim.get("status") == "ok"
    for i, v in enumerate(vars_plot):
        ax = axs[i]
        sub = sub_all[sub_all["variable"] == v]
        if "sigma" in sub.columns and sub["sigma"].notna().any():
            ax.errorbar(sub["t_h"], sub["value"], yerr=sub["sigma"],
                        fmt="o", color="tab:red", ms=4, capsize=2, label="exp")
        else:
            ax.plot(sub["t_h"], sub["value"], "o", color="tab:red", ms=4, label="exp")
        if ok:
            y_sim = sim["obs"].get(v)
            if y_sim is not None and np.any(np.isfinite(y_sim)):
                ax.plot(sim["t_h"], y_sim, "-", color="tab:gray", lw=1.8, label="sim nominal")
        ax.set_title(v, fontsize=9)
        ax.set_xlabel("t (h)", fontsize=8)
        ax.grid(alpha=0.3)
        if i == 0:
            ax.legend(loc="best", fontsize=7)
    for j in range(len(vars_plot), len(axs)):
        axs[j].set_visible(False)

    meta = experiments_meta[dataset_id]
    ic = meta["initial_conditions"]
    Tmin = meta["T_profile"]["T_K"].min() - 273.15 if len(meta["T_profile"]) else float("nan")
    Tmax = meta["T_profile"]["T_K"].max() - 273.15 if len(meta["T_profile"]) else float("nan")
    fig.suptitle(
        f"[{dataset_id}] nominal vs exp | t_end={meta['t_end']:.0f}h | "
        f"T={Tmin:.1f}-{Tmax:.1f}C | pulsos={len(meta['pulse_schedule'])} | "
        f"IC X={ic['X_gDW_L']:.3f} N={ic['N_gN_L']:.4f}",
        fontsize=9, y=1.02)
    plt.tight_layout()
    plt.show()
    return sim


def pre_pso_verification(datasets=None, theta_unit=None):
    """Loop sobre datasets: simular nominal con IC reales + plot + RMSE."""
    if experiments_meta is None or not experiments_meta:
        print("[VERIF] experiments_meta vacio (sin datos reales). Skip.")
        return {}
    if datasets is None:
        datasets = list(experiments_meta.keys())
    if theta_unit is None:
        theta_unit = theta_nominal_unit

    print("=" * 80)
    print(f"VERIFICACION PRE-PSO: theta_nominal vs datos ({len(datasets)} ensayos)")
    print("IC + T-profile + pulsos ruteados al solver. Fase 5 (multi-exp objective) abajo.")
    print("=" * 80)

    results = {}
    rmse_rows = []
    for ds in datasets:
        print(f"\n[VERIF] {ds} ...")
        sim = forward_simulate_experiment(theta_unit, ds)
        results[ds] = sim
        if sim.get("status") != "ok":
            print(f"  FAIL: {sim.get('error', '?')}")
            continue
        rmse_map = _residuals_for_dataset(sim, ds) or {}
        row = {"dataset_id": ds}
        for v, r in rmse_map.items():
            row[v] = r
        rmse_rows.append(row)
        # Plot
        plot_nominal_vs_experiment(ds, sim=sim)

    if rmse_rows:
        rmse_df = pd.DataFrame(rmse_rows).set_index("dataset_id")
        print("\n[VERIF] RMSE_norm por ensayo y observable (baseline con nominal):")
        print(rmse_df.round(3).to_string())
    return results


# --- Ejecucion -------------------------------------------------------------------------
if experiments_meta and ("forward_simulate" in globals()):
    _verif_sims = pre_pso_verification()
else:
    print("[VERIF] Pipeline incompleto (sin experiments_meta o sin forward_simulate). "
          "Ejecuta celdas 61 y 65 primero.")
