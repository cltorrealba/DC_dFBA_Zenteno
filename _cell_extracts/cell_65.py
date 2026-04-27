# ======================================================================================
# MODULO 3 — Forward simulation wrapper
# ======================================================================================
import warnings
import traceback

# --- Configuracion del forward solver --------------------------------------------------
NFE_OPT   = 12   # nfe para fase de optimizacion (reducir para acelerar PSO)
NCP_OPT   = 3
TEND_OPT  = 72.0
NITER_OPT = 2

# --- Toggles de verbosidad del forward solver -----------------------------------------
# OPT_VERBOSE_SOLVER = False  -> suprime stdout, plots Y tablas del solver dFBA en CADA
#                                 evaluacion (recomendado para PSO).
# OPT_VERBOSE_SOLVER = True   -> deja pasar todo (modo debug manual, 1 llamada).
OPT_VERBOSE_SOLVER = False

# Alias retro-compatible (lo lee forward_simulate si silent=None).
FORWARD_SILENT = not OPT_VERBOSE_SOLVER

# Mapeo: observable -> columna en states_df (o proceso especial)
OBS_TO_COL = {
    "X_gDW_L":             ("state",  "X_gDW_L"),
    "G_g_L":               ("state",  "G_g_L"),
    "F_g_L":               ("state",  "F_g_L"),
    "E_g_L":               ("state",  "E_g_L"),
    "N_gN_L":              ("state",  "N_gN_L"),
    # Acetate esters (integracion via flujos de colocacion)
    "isoamyl_acetate_mg_L":        ("aroma",  "r_1862"),   # isoamyl acetate
    "ethyl_acetate_mg_L":  ("aroma",  "r_1765"),   # ethyl acetate
}

# Mapa MW para integracion de aromas (g/mmol).
AROMA_MW_MAP_OPT = {
    "r_1862": 0.130186,   # isoamyl acetate
    "r_1765": 0.088106,   # ethyl acetate
    "r_1867": 0.116160,   # isobutyl acetate
    "r_2000": 0.164200,   # phenethyl acetate
}

# RIDs de aromas extra que el solver balanced debe extraer explicitamente
# (no estan en found_aromas porque ese dict mapea los ALCOHOLES, no los acetatos).
EXTRA_AROMA_RIDS = ["r_1862", "r_1765", "r_1867", "r_2000"]

# Mapeo rid -> nombre observable (para aroma IC desde AROMA_IC_MG_L).
AROMA_RID_TO_OBS = {
    "r_1862": "isoamyl_acetate_mg_L",
    "r_1765": "ethyl_acetate_mg_L",
}

# Inspecciona nombres sospechosos
suspect = [c.name for c in model.constraints if 'phi' in c.name or 'ester' in c.name.lower() or 'acetate' in c.name.lower()]
print(suspect)


def _integrate_aroma_mg_L(states_df, flux_df, aroma_rid, mw_g_per_mmol):
    """
    Integra flujo de aroma (mmol/gDW/h) a concentracion mg/L por trapezoidal,
    replicando la logica de celda 46:
        dC/dt = 1000 * MW * v * X       =>   C(t_{k+1}) = C(t_k) + MW*1000 * v_k * X_k * dt
    Usa el valor del flujo en el punto de colocacion anterior (explicit Euler,
    consistente con celda 46). Para mayor precision se puede cambiar a trapezoidal.
    """
    # Primario: columna directa v_aroma_<rid>
    flux_col = None
    cand_rid = f"v_aroma_{aroma_rid}"
    if cand_rid in flux_df.columns:
        flux_col = cand_rid
    # Secundario: invertir found_aromas (key -> rid) y probar v_aroma_<key>
    if flux_col is None:
        _fa = globals().get("found_aromas", {}) or {}
        for _ak, _rid in _fa.items():
            if _rid == aroma_rid:
                cand_ak = f"v_aroma_{_ak}"
                if cand_ak in flux_df.columns:
                    flux_col = cand_ak
                    break
    # Terciario: compat con nomenclatura plana v_<key> (solvers antiguos)
    if flux_col is None:
        _fa = globals().get("found_aromas", {}) or {}
        for _ak, _rid in _fa.items():
            if _rid == aroma_rid:
                cand_bare = f"v_{_ak}"
                if cand_bare in flux_df.columns:
                    flux_col = cand_bare
                    break
    if flux_col is None:
        return np.full(len(states_df), np.nan)

    t_state = states_df["t_h"].to_numpy(dtype=float)
    X_state = states_df["X_gDW_L"].to_numpy(dtype=float)
    t_flux  = flux_df["t_h"].to_numpy(dtype=float)
    v_flux  = np.maximum(0.0, flux_df[flux_col].to_numpy(dtype=float))

    C = np.zeros(len(t_state))
    # IC de aroma desde AROMA_IC_MG_L (seteado por set_process_inputs si hay datos).
    _obs_name = AROMA_RID_TO_OBS.get(aroma_rid)
    if _obs_name is not None:
        _ic_map = globals().get("AROMA_IC_MG_L", {}) or {}
        C[0] = float(_ic_map.get(_obs_name, 0.0))
    for k in range(len(t_state) - 1):
        dt = t_state[k + 1] - t_state[k]
        # v en el punto de colocacion asociado (el k-esimo de flux_df suele alinear con k+1 del estado)
        vk = float(v_flux[k]) if k < len(v_flux) else 0.0
        Xk = float(X_state[k])
        C[k + 1] = C[k] + 1000.0 * mw_g_per_mmol * vk * Xk * dt
    return C


# Suprimir side-effects de plotting + stdout + rich-display de _run_dfba_colloc*
import matplotlib
import matplotlib.pyplot as plt
from contextlib import contextmanager

try:
    from IPython.utils.capture import capture_output as _ipy_capture
    _HAS_IPY = True
except ImportError:
    _HAS_IPY = False


@contextmanager
def _silence_all():
    """
    Silencia completamente el output del solver:
      - plt.show() => no renderiza
      - cualquier figura creada dentro del bloque => cerrada
      - stdout / stderr del solver => capturado (stringbuffer desechado)
      - display() / rich-representacion (tablas Jupyter) => capturado via IPython
    """
    _orig_show = plt.show
    plt.show = lambda *a, **kw: None
    figs_before = set(plt.get_fignums())

    if _HAS_IPY:
        with _ipy_capture(stdout=True, stderr=True, display=True):
            try:
                yield
            finally:
                figs_after = set(plt.get_fignums())
                for fnum in (figs_after - figs_before):
                    plt.close(fnum)
                plt.show = _orig_show
    else:
        import io, contextlib
        buf_out, buf_err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(buf_out), contextlib.redirect_stderr(buf_err):
            try:
                yield
            finally:
                figs_after = set(plt.get_fignums())
                for fnum in (figs_after - figs_before):
                    plt.close(fnum)
                plt.show = _orig_show


# Alias retro-compatible
_suppress_plots = _silence_all


def forward_simulate(theta_unit, *, nfe=NFE_OPT, ncp=NCP_OPT, t_end=TEND_OPT,
                    n_iter=NITER_OPT, verbose=False, silent=None):
    """
    Corre un forward dFBA con los parametros theta_unit (espacio [0,1]^n).
    Retorna dict con trayectorias de los observables.

    silent=None (default): usa FORWARD_SILENT global.
    silent=True            : suprime plots y prints del solver dFBA (critico para PSO).
    silent=False           : deja pasar plots y prints (debug).
    """
    if silent is None:
        silent = bool(FORWARD_SILENT)
    phys = apply_theta(theta_unit, verbose=verbose)

    import contextlib
    _ctx = _silence_all() if silent else contextlib.nullcontext()

    try:
     with _ctx:
                # Preferimos el solver balanceado (celda 46); fallback al basico si no existe.
                if "_run_dfba_colloc_balanced" in globals():
                    res = _run_dfba_colloc_balanced(
                        scenario_name="opt_forward",
                        aerobic_mode=False,
                        o2_init_g_l=0.0,
                        nfe=nfe, ncp=ncp, t_end=t_end, n_iter=n_iter,
                        growth_lock_frac=1.0, w_atpm_growth=0.0, w_aa_growth=0.0,
                        use_pfba_mix=True,
                    )
                elif "_run_dfba_colloc" in globals():
                    res = _run_dfba_colloc(
                        scenario_name="opt_forward",
                        aerobic_mode=False, o2_init_g_l=0.0,
                        nfe=nfe, ncp=ncp, t_end=t_end, n_iter=n_iter,
                    )
                else:
                    raise RuntimeError("No hay solver dFBA disponible (_run_dfba_colloc*).")
    except Exception as e:
        # Log de los primeros N fallos para diagnostico (evita flood en PSO).
        _n = globals().get("_FORWARD_FAIL_LOG_COUNT", 0)
        _lim = globals().get("_FORWARD_FAIL_LOG_LIMIT", 5)
        if _n < _lim:
            import traceback as _tb
            try:
                import sys as _sys
                _sys.__stdout__.write(
                    f"[forward#FAIL {_n+1}/{_lim}] {type(e).__name__}: {e}\n")
                _sys.__stdout__.flush()
                if verbose:
                    _sys.__stdout__.write(_tb.format_exc() + "\n")
                    _sys.__stdout__.flush()
            except Exception:
                print(f"[forward#FAIL {_n+1}/{_lim}] {type(e).__name__}: {e}")
            globals()["_FORWARD_FAIL_LOG_COUNT"] = _n + 1
        return {"status": "fail", "error": str(e), "theta_phys": phys}

    states_df = res.get("states_df")
    flux_df   = res.get("flux_df")
    if states_df is None or len(states_df) == 0:
        return {"status": "fail", "error": "states_df vacio", "theta_phys": phys}

    # --- Construir diccionario de observables ----------------------------------------
    t_h = states_df["t_h"].to_numpy(dtype=float)
    obs = {}

    for var in TARGET_OBSERVABLES:
        kind, key = OBS_TO_COL.get(var, (None, None))
        if kind == "state":
            obs[var] = (states_df[key].to_numpy(dtype=float)
                        if key in states_df.columns else np.full(len(t_h), np.nan))
        elif kind == "aroma":
            mw = AROMA_MW_MAP_OPT.get(key, 0.1)
            obs[var] = _integrate_aroma_mg_L(states_df, flux_df, key, mw)
        else:
            obs[var] = np.full(len(t_h), np.nan)

    return {
        "status":     "ok",
        "t_h":        t_h,
        "obs":        obs,
        "states_df":  states_df,
        "flux_df":    flux_df,
        "theta_phys": phys,
    }


# --- Smoke test con parametros nominales -----------------------------------------------
print("[FORWARD] Smoke test con theta nominal (nfe={}, ncp={})...".format(NFE_OPT, NCP_OPT))
sim_nom = forward_simulate(theta_nominal_unit, verbose=False)
print(f"[FORWARD] status = {sim_nom['status']}")
if sim_nom["status"] == "ok":
    for v, arr in sim_nom["obs"].items():
        finite = arr[np.isfinite(arr)]
        if len(finite):
            print(f"  {v:18s} t_final={arr[-1]:10.4f}  min={finite.min():.4f}  max={finite.max():.4f}")
        else:
            print(f"  {v:18s} (sin datos)")


# ======================================================================================
# MODULO 3b — Instalacion de restricciones soft de esteres (persistentes sobre model)
# ======================================================================================
# Estas constraints se aplican DIRECTAMENTE al `model` global (no al snapshot). Como
# el solver _fba_at_state_colloc_balanced usa `with model as mtmp`, las constraints
# permanecen presentes en cada LP del forward.
#
# Se actualizan antes de cada forward_simulate via el hook _sync_ester_soft_constraints()
# registrado en apply_theta.

ESTER_CONS_NAMES = []   # nombres de constraints registradas (para poder removerlas/renovar)


def _remove_ester_soft_constraints(sweep_by_pattern=True):
    """
    Remueve constraints soft de esteres en AMBAS capas (cobra + optlang).

    Estrategia:
      1. Junta nombres candidatos desde `ESTER_CONS_NAMES`, `model.constraints` y
         `model.solver.constraints` (optlang).
      2. Filtra por patron si sweep_by_pattern=True.
      3. Remueve via `model.solver.remove(name)` (optlang) — capa autoritativa.
         Si falla, intenta `model.remove_cons_vars` como fallback.

    Patrones barridos:
      - "pair_<ester>__<alcohol>_phi<val>"   (celda 41, pairwise)
      - "ester_coupling_phi_<val>"           (celda 40, global)
      - "ethyl_acetate_soft_phi<val>"        (celda 42, EA soft)
    """
    import re as _re
    global ESTER_CONS_NAMES

    patterns = [
        _re.compile(r"^pair_.+__.+_phi"),
        _re.compile(r"^ester_coupling_phi_"),
        _re.compile(r"^ethyl_acetate_soft_phi"),
    ]

    # 1) Reunir nombres candidatos desde todas las fuentes
    candidates = set(ESTER_CONS_NAMES)

    try:
        candidates.update([c.name for c in model.constraints])
    except Exception:
        pass

    try:
        # optlang container (autoritativo)
        candidates.update(list(model.solver.constraints.keys()))
    except Exception:
        pass

    # 2) Filtrar por patron
    if sweep_by_pattern:
        targets = [n for n in candidates if any(p.match(n) for p in patterns)]
    else:
        targets = [n for n in candidates if n in ESTER_CONS_NAMES]

    # 3) Remover en capa optlang (autoritativa)
    removed = 0
    for name in targets:
        # Intento 1: solver (optlang)
        try:
            if name in model.solver.constraints:
                model.solver.remove(name)
                removed += 1
                continue
        except Exception:
            pass
        # Intento 2: cobra wrapper
        try:
            if name in model.constraints:
                model.remove_cons_vars([model.constraints[name]])
                removed += 1
        except Exception:
            pass

    ESTER_CONS_NAMES = []
    return removed


def _install_ester_soft_constraints(phi_pair_val=None, phi_ea_val=None, verbose=False):
    """
    Instala las restricciones soft de esteres sobre `model` usando los helpers
    existentes en celdas 41/42:
      - _add_pairwise_ester_coupling_constraints(model, df)    => ester/alcohol por par
      - _add_ethyl_acetate_soft_constraint(model, phi_ea)      => ethyl acetate / ethanol

    phi_pair_val: valor actual de PAIRWISE_PHI_STATIC. Si None => lee global.
    phi_ea_val  : valor actual de PHI_ETHYL_ACETATE_STATIC. Si None => lee global.
    """
    global ESTER_CONS_NAMES

    # 1) Limpiar cualquier instalacion previa
    _remove_ester_soft_constraints()

    phi_pair_val = float(phi_pair_val if phi_pair_val is not None
                         else globals().get("PAIRWISE_PHI_STATIC", 0.0))
    phi_ea_val   = float(phi_ea_val   if phi_ea_val   is not None
                         else globals().get("PHI_ETHYL_ACETATE_STATIC", 0.0))

    # 2) Pairwise ester/alcohol (requiere ester_alcohol_pairwise_df de celda 41)
    pw_df = globals().get("ester_alcohol_pairwise_df", None)
    fn_pair = globals().get("_add_pairwise_ester_coupling_constraints", None)
    if pw_df is not None and fn_pair is not None and phi_pair_val > 0:
        pw_work = pw_df.copy()
        if "pair_active" in pw_work.columns:
            pw_work.loc[pw_work["pair_active"], "phi_pair"] = phi_pair_val
        try:
            cons_added, _am, _pm = fn_pair(model, pw_work)
            for c in (cons_added or []):
                ESTER_CONS_NAMES.append(c.name)
        except Exception as e:
            if verbose:
                print(f"[esters] fallo instalar pairwise: {e}")

    # 3) Ethyl acetate soft (requiere _add_ethyl_acetate_soft_constraint de celda 42)
    fn_ea = globals().get("_add_ethyl_acetate_soft_constraint", None)
    # Asegurar que INCLUDE_ETHYL_ACETATE_SOFT_COUPLING=True globalmente
    globals()["INCLUDE_ETHYL_ACETATE_SOFT_COUPLING"] = True
    if fn_ea is not None and phi_ea_val > 0:
        try:
            ea_cons, ea_added = fn_ea(model, phi_ea_val)
            if ea_added and ea_cons is not None:
                ESTER_CONS_NAMES.append(ea_cons.name)
        except Exception as e:
            if verbose:
                print(f"[esters] fallo instalar ethyl acetate soft: {e}")

    if verbose:
        print(f"[esters] instaladas {len(ESTER_CONS_NAMES)} constraints "
              f"(phi_pair={phi_pair_val:.4g}, phi_ea={phi_ea_val:.4g}).")


def _sync_ester_soft_constraints():
    """Hook llamado por apply_theta(). Reinstala con los phi actuales."""
    _install_ester_soft_constraints(verbose=False)


# Instalacion inicial con valores actuales (solo si los helpers estan disponibles)
if ("_add_pairwise_ester_coupling_constraints" in globals()
    or "_add_ethyl_acetate_soft_constraint" in globals()):
    _install_ester_soft_constraints(verbose=True)
else:
    print("[esters] Helpers de constraints no encontrados. Ejecuta celdas 41 y 42 antes "
          "de correr el PSO para habilitar calibracion de esteres.")


# ======================================================================================
# Diagnostico rapido del pipeline
# ======================================================================================
def diagnose_forward():
    """Imprime checklist del estado previo a correr PSO."""
    import sys as _sys
    checks = []
    checks.append(("model presente", "model" in globals()))
    checks.append(("_run_dfba_colloc_balanced disponible", "_run_dfba_colloc_balanced" in globals()))
    checks.append(("_run_dfba_colloc disponible", "_run_dfba_colloc" in globals()))
    checks.append(("_add_pairwise_ester_coupling_constraints disponible",
                   "_add_pairwise_ester_coupling_constraints" in globals()))
    checks.append(("_add_ethyl_acetate_soft_constraint disponible",
                   "_add_ethyl_acetate_soft_constraint" in globals()))
    checks.append(("ester_alcohol_pairwise_df presente",
                   "ester_alcohol_pairwise_df" in globals()))
    checks.append(("states_df_phase presente", "states_df_phase" in globals()))
    checks.append(("exp_df presente", "exp_df" in globals()))
    checks.append(("theta_nominal_unit definido", "theta_nominal_unit" in globals()))

    ok_all = True
    print("=" * 70)
    print("DIAGNOSTICO PIPELINE PSO")
    print("=" * 70)
    for label, ok in checks:
        mark = "OK" if ok else "FALTA"
        if not ok:
            ok_all = False
        print(f"  [{mark}] {label}")

    # Intentar corrida nominal con verbose de error
    print("\n--- Smoke test forward con theta_nominal (silent=False) ---")
    globals()["_FORWARD_FAIL_LOG_COUNT"] = 0  # reset contador
    sim = forward_simulate(theta_nominal_unit, silent=True, verbose=True)
    print(f"  status = {sim.get('status')}")
    if sim.get("status") != "ok":
        print(f"  error  = {sim.get('error', '?')}")
    else:
        for k, arr in sim["obs"].items():
            import numpy as _np
            finite = arr[_np.isfinite(arr)] if arr is not None else []
            print(f"  {k:24s} final={arr[-1] if len(arr) else _np.nan:10.4g}  "
                  f"min={finite.min() if len(finite) else _np.nan:10.4g}  "
                  f"max={finite.max() if len(finite) else _np.nan:10.4g}")

    if not ok_all:
        print("\n[!!] Hay dependencias faltantes. Revisa las celdas 41/42/44/45/46/61/63/65.")
    elif sim.get("status") != "ok":
        print("\n[!!] forward_simulate falla con nominal. Revisa el traceback arriba.")
    else:
        print("\n[OK] Pipeline listo para PSO.")
    return sim


# ======================================================================================
# MODULO 3c - Generacion de datos sinteticos desde forward_simulate(theta_truth)
# ======================================================================================
# Si exp_df fue diferido en celda 61 (no hay CSV real), aqui lo generamos usando la
# simulacion nominal: theta_truth = theta_nominal (incluye phi_pair=0.08, phi_ea=1e-4).
# Esto garantiza un "ground truth" coherente con el modelo cell-43-style (esters > 0).

def _generate_synthetic_exp_df(sim_truth, n_points=None, noise_frac=None,
                                noise_floor=None, seed=None, dataset_id="SYNTH_TRUTH"):
    """Convierte una simulacion forward en un DataFrame long-form de pseudo-datos."""
    n_points    = int(n_points    if n_points    is not None else SYNTH_N_POINTS)
    noise_frac  = float(noise_frac  if noise_frac  is not None else SYNTH_NOISE_FRAC)
    noise_floor = float(noise_floor if noise_floor is not None else SYNTH_NOISE_FLOOR)
    seed        = int(seed        if seed        is not None else SYNTH_SEED)

    if sim_truth.get("status") != "ok":
        raise RuntimeError(f"sim_truth fallo: {sim_truth.get('error', '?')}")

    t_sim = sim_truth["t_h"]
    obs   = sim_truth["obs"]
    t_exp = np.linspace(float(t_sim[0]), float(t_sim[-1]), n_points)
    rng   = np.random.default_rng(seed)

    rows = []
    for v in TARGET_OBSERVABLES:
        arr = obs.get(v)
        if arr is None or not np.any(np.isfinite(arr)):
            print(f"  [SYNTH] {v}: sin datos en sim_truth, skip.")
            continue
        mask = np.isfinite(arr)
        y_clean = np.interp(t_exp, t_sim[mask], arr[mask])
        for t_i, y_i in zip(t_exp, y_clean):
            sigma_i = max(noise_frac * abs(y_i), noise_floor)
            noise   = rng.normal(0.0, sigma_i)
            rows.append({
                "dataset_id": dataset_id,
                "t_h":        float(t_i),
                "variable":   v,
                "value":      float(y_i + noise),
                "sigma":      float(sigma_i),
            })
    return pd.DataFrame(rows)


# Decidir si debemos generar sinteticos
_need_synth = (
    ("exp_df" not in globals())
    or (exp_df is None)
    or (len(exp_df) == 0)
    or globals().get("FORCE_SYNTH_REGEN", False)
)

if _need_synth and globals().get("USE_SYNTHETIC_FROM_FORWARD", True):
    print("\n" + "=" * 80)
    print("[SYNTH] Generando exp_df sintetico desde forward_simulate(theta_nominal)")
    print(f"[SYNTH] phi_pair_truth = {globals().get('PAIRWISE_PHI_STATIC', 0.08)}, "
          f"phi_ea_truth = {globals().get('PHI_ETHYL_ACETATE_STATIC', 1e-4)}")
    print("=" * 80)
    if sim_nom.get("status") != "ok":
        raise RuntimeError("sim_nom no OK; no se puede generar sintetico. "
                           "Revisa forward_simulate / instalacion de constraints.")
    exp_df = _generate_synthetic_exp_df(sim_nom)
    print(f"[SYNTH] Generadas {len(exp_df)} obs sobre {exp_df['variable'].nunique()} variables.")
    print("[SYNTH] Resumen (n por variable):")
    print(exp_df.groupby("variable").size().to_string())
    print(f"[SYNTH] Rango temporal: {exp_df['t_h'].min():.1f} - {exp_df['t_h'].max():.1f} h")
    # Sanity: aromas > 0?
    for v in ["isoamyl_acetate_mg_L", "ethyl_acetate_mg_L"]:
        sub = exp_df[exp_df["variable"] == v]["value"]
        if len(sub):
            print(f"[SYNTH] {v}: min={sub.min():.4g}, max={sub.max():.4g}, mean={sub.mean():.4g}")
