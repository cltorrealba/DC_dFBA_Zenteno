# ======================================================================================
# MODULO 5 - Problema pymoo + PSO (con checkpointing + Ctrl+C)
# ======================================================================================
import sys
import time
import json as _json
from pathlib import Path as _Path

try:
    from pymoo.core.problem import ElementwiseProblem
    from pymoo.core.callback import Callback
    from pymoo.algorithms.soo.nonconvex.pso import PSO
    from pymoo.optimize import minimize
    from pymoo.termination import get_termination
    _PYMOO_OK = True
except ImportError as _e:
    _PYMOO_OK = False
    print(f"[PSO] pymoo no disponible: {_e}")


# --- Configuracion -------------------------------------------------------------------
PSO_POP_SIZE = 30
PSO_N_GEN    = 50
PSO_SEED     = 42
PSO_VERBOSE  = True
PSO_PROGRESS_EVERY = 0   # 0 = no imprimir per-eval (deja que pymoo muestre por gen)

# Checkpointing
PSO_CHECKPOINT_DIR  = _Path("out") / "pso_checkpoints"
PSO_CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
PSO_CHECKPOINT_BEST = PSO_CHECKPOINT_DIR / "pso_best_latest.json"
PSO_CHECKPOINT_GEN  = PSO_CHECKPOINT_DIR / "pso_state_latest.json"


# Estado global compartido (sobrevive a Ctrl+C / crash)
_pso_state = {
    "best_F":    float("inf"),
    "best_X":    None,        # vector unit (lista)
    "best_phys": None,        # dict nombre->valor fisico
    "n_eval":    0,
    "n_gen":     0,
    "elapsed_s": 0.0,
    "history":   [],          # uno por evaluacion: {eval, J, dt}
    "gen_log":   [],          # uno por generacion: {gen, n_eval, F_min, F_avg}
}


def _save_checkpoint_best():
    """Guarda solo el best-so-far (escritura barata, llamada en cada nuevo best)."""
    payload = {
        "best_F":    _pso_state["best_F"],
        "best_X":    _pso_state["best_X"],
        "best_phys": _pso_state["best_phys"],
        "n_eval":    _pso_state["n_eval"],
        "n_gen":     _pso_state["n_gen"],
    }
    try:
        with PSO_CHECKPOINT_BEST.open("w", encoding="utf-8") as f:
            _json.dump(payload, f, indent=2, default=float)
    except Exception as e:
        sys.__stdout__.write(f"[CKPT] WARN: no pude guardar best ({e})\n")


def _save_checkpoint_full():
    """Guarda estado completo (history + gen_log). Mas pesado, llamada por gen."""
    try:
        with PSO_CHECKPOINT_GEN.open("w", encoding="utf-8") as f:
            _json.dump(_pso_state, f, indent=2, default=float)
    except Exception as e:
        sys.__stdout__.write(f"[CKPT] WARN: no pude guardar state ({e})\n")


class ParamEstProblem(ElementwiseProblem):
    """Problema 1-objetivo, caja unitaria [0,1]^n. Actualiza _pso_state."""
    def __init__(self, n_var, **kwargs):
        super().__init__(n_var=n_var, n_obj=1, n_constr=0,
                         xl=np.zeros(n_var), xu=np.ones(n_var), **kwargs)
        self._t_start = time.time()

    def _evaluate(self, x, out, *args, **kwargs):
        t0 = time.time()
        x_arr = np.asarray(x, dtype=float)
        J = objective(x_arr)
        dt = time.time() - t0
        out["F"] = J

        _pso_state["n_eval"]    += 1
        _pso_state["elapsed_s"]  = time.time() - self._t_start
        _pso_state["history"].append({"eval": _pso_state["n_eval"], "J": float(J), "dt_s": dt})

        if np.isfinite(J) and J < _pso_state["best_F"]:
            _pso_state["best_F"]    = float(J)
            _pso_state["best_X"]    = x_arr.tolist()
            try:
                _pso_state["best_phys"] = {k: float(v) for k, v in
                                            theta_to_physical(x_arr).items()}
            except Exception:
                _pso_state["best_phys"] = None
            _save_checkpoint_best()

        if PSO_PROGRESS_EVERY > 0 and (_pso_state["n_eval"] % PSO_PROGRESS_EVERY == 0):
            try:
                sys.__stdout__.write(
                    f"[PSO] eval {_pso_state['n_eval']:4d}  J={J:12.5g}  "
                    f"J_best={_pso_state['best_F']:12.5g}  dt={dt:5.1f}s\n")
                sys.__stdout__.flush()
            except Exception:
                pass


class _CheckpointCallback(Callback):
    """Callback al final de cada generacion: log + dump full state."""
    def notify(self, algorithm):
        _pso_state["n_gen"] += 1
        try:
            F = np.asarray(algorithm.pop.get("F"), dtype=float).ravel()
            F_min = float(np.nanmin(F))
            F_avg = float(np.nanmean(F))
        except Exception:
            F_min, F_avg = float("nan"), float("nan")
        _pso_state["gen_log"].append({
            "gen":   _pso_state["n_gen"],
            "n_eval": _pso_state["n_eval"],
            "F_min": F_min,
            "F_avg": F_avg,
            "best_so_far": _pso_state["best_F"],
        })
        _save_checkpoint_full()


def _print_best_summary():
    """Imprime el theta best-so-far guardado en _pso_state."""
    if _pso_state["best_X"] is None:
        print("[PSO] (sin best registrado)")
        return
    print(f"\n[PSO] best-so-far: J = {_pso_state['best_F']:.6g}  "
          f"(eval={_pso_state['n_eval']}, gen={_pso_state['n_gen']})")
    print(f"[PSO] checkpoint best   : {PSO_CHECKPOINT_BEST}")
    print(f"[PSO] checkpoint state  : {PSO_CHECKPOINT_GEN}")
    if _pso_state["best_phys"]:
        print("[PSO] theta_best (fisico):")
        for name, val in _pso_state["best_phys"].items():
            try:
                spec = next(p for p in active_specs() if p.name == name)
                flag = ""
                rng = max(abs(spec.ub - spec.lb), 1e-12)
                if abs(val - spec.lb) / rng < 1e-3:
                    flag = "  <-- LB"
                elif abs(val - spec.ub) / rng < 1e-3:
                    flag = "  <-- UB"
                print(f"  {name:32s} = {val:12.6g}   (nom={spec.nominal:.4g}){flag}")
            except StopIteration:
                print(f"  {name:32s} = {val:12.6g}")


def reset_pso_state():
    """Limpia el estado global (llama antes de cada corrida nueva)."""
    _pso_state.update({
        "best_F": float("inf"), "best_X": None, "best_phys": None,
        "n_eval": 0, "n_gen": 0, "elapsed_s": 0.0,
        "history": [], "gen_log": [],
    })


def run_pso(pop_size=PSO_POP_SIZE, n_gen=PSO_N_GEN, seed=PSO_SEED,
            verbose=PSO_VERBOSE, skip_sanity=False, resume=False):
    """
    Corre PSO con checkpoint + Ctrl+C-safe.
      - Estado se persiste en _pso_state (memoria) y en out/pso_checkpoints/*.json (disco).
      - Ctrl+C: captura KeyboardInterrupt, retorna best-so-far en lugar de propagar.
      - resume=True: NO resetea _pso_state (uso para continuar tras Ctrl+C, aunque
        pymoo internamente arranca de cero; se preserva el best_X historico).
    """
    if not _PYMOO_OK:
        raise RuntimeError("pymoo no esta instalado.")

    if not skip_sanity:
        print("[PSO] Sanity check con theta_nominal ...")
        globals()["_FORWARD_FAIL_LOG_COUNT"] = 0
        sim_sanity = forward_simulate(theta_nominal_unit, silent=True)
        if sim_sanity.get("status") != "ok":
            print(f"[PSO] ABORT: forward_simulate falla con nominal: "
                  f"{sim_sanity.get('error', '?')}")
            raise RuntimeError("Forward falla con nominal.")
        _, J_sanity = residuals_by_observable(sim_sanity)
        print(f"[PSO] OK: J(nominal) = {J_sanity:.6g}")

    if not resume:
        reset_pso_state()

    problem     = ParamEstProblem(n_var=N_PARAMS_ACTIVE)
    algorithm   = PSO(pop_size=pop_size)
    termination = get_termination("n_gen", n_gen)
    callback    = _CheckpointCallback()

    print(f"[PSO] Iniciando: n_var={N_PARAMS_ACTIVE}, pop={pop_size}, n_gen={n_gen}")
    print(f"[PSO] Checkpoints en {PSO_CHECKPOINT_DIR}")
    print("[PSO] Tip: Ctrl+C detiene y conserva el best-so-far.")

    globals()["_FORWARD_FAIL_LOG_COUNT"] = 0
    t_start = time.time()
    res = None
    interrupted = False
    try:
        res = minimize(
            problem, algorithm, termination,
            seed=seed, verbose=verbose, save_history=False,
            callback=callback,
        )
    except KeyboardInterrupt:
        interrupted = True
        print("\n[PSO] *** Ctrl+C detectado: deteniendo PSO y conservando best-so-far ***")
    except Exception as e:
        print(f"\n[PSO] *** Error inesperado: {type(e).__name__}: {e} ***")
        print("[PSO] Conservando best-so-far en _pso_state / checkpoint.")

    elapsed = time.time() - t_start
    print(f"[PSO] Terminado en {elapsed:.1f} s, {_pso_state['n_eval']} evaluaciones, "
          f"{_pso_state['n_gen']} generaciones.")

    _print_best_summary()

    # Construir resultado con el best-so-far (no con res.X de pymoo, que puede
    # ser None si Ctrl+C interrumpio antes de la primera generacion completa).
    theta_star_unit = (np.asarray(_pso_state["best_X"], dtype=float)
                       if _pso_state["best_X"] is not None else None)
    theta_star_phys = (theta_to_physical(theta_star_unit)
                       if theta_star_unit is not None else None)

    return {
        "theta_star_unit": theta_star_unit,
        "theta_star_phys": theta_star_phys,
        "J_star":          _pso_state["best_F"],
        "pymoo_result":    res,
        "problem":         problem,
        "elapsed_s":       elapsed,
        "interrupted":     interrupted,
        "n_eval":          _pso_state["n_eval"],
        "n_gen":           _pso_state["n_gen"],
        "checkpoint_best": str(PSO_CHECKPOINT_BEST),
        "checkpoint_full": str(PSO_CHECKPOINT_GEN),
    }


def load_pso_checkpoint(best_only=True):
    """Recupera el best-so-far desde disco (utilidad post-crash)."""
    p = PSO_CHECKPOINT_BEST if best_only else PSO_CHECKPOINT_GEN
    if not p.exists():
        print(f"[CKPT] No existe {p}")
        return None
    with p.open("r", encoding="utf-8") as f:
        data = _json.load(f)
    print(f"[CKPT] Cargado {p}: best_F={data.get('best_F')}, "
          f"n_eval={data.get('n_eval')}, n_gen={data.get('n_gen')}")
    return data


# --- Ejecucion -----------------------------------------------------------------------
pso_result = run_pso(pop_size=PSO_POP_SIZE, n_gen=PSO_N_GEN)

print(f"\n[PSO] Dimension: {N_PARAMS_ACTIVE} | Pop: {PSO_POP_SIZE} | Gen: {PSO_N_GEN}")
print(f"[PSO] LPs estimados por eval: ~{NFE_OPT * NCP_OPT * NITER_OPT}")
