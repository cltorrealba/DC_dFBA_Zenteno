# ======================================================================================
# MODULO 2 — Definicion del vector de parametros
# ======================================================================================
from dataclasses import dataclass, field
from typing import List, Dict


@dataclass
class ParamSpec:
    name: str            # nombre de la global a mutar
    nominal: float       # valor nominal
    lb: float            # bound inferior (en espacio fisico)
    ub: float            # bound superior (en espacio fisico)
    log_scale: bool = False   # optimizar en log10 (True para K*, phi)
    group: str = "kinetic"    # kinetic / aroma / turnover / init
    enabled: bool = True      # si False, se mantiene nominal y no entra al vector


# --- Especificacion maestra -----------------------------------------------------------
# Nota: los "nominal" intentan leerse desde las globales actuales, con fallback al paper.

def _g(name, default):
    return float(globals().get(name, default))


PARAM_SPECS: List[ParamSpec] = [
    # === Grupo 1: Cinetica Zenteno ===
    ParamSpec("MU0_nom",     _g("MU0_nom",     0.18),   0.08,   0.35,  False, "kinetic", True),
    ParamSpec("YXN_nom",     _g("YXN_nom",    19.69),   5.0,   40.0,   False, "kinetic", True),
    ParamSpec("YXG_nom",     _g("YXG_nom",     1.60),   0.3,    3.0,   False, "kinetic", True),
    ParamSpec("YXF_nom",     _g("YXF_nom",     1.60),   0.3,    3.0,   False, "kinetic", True),
    ParamSpec("YEG_nom",     _g("YEG_nom",     0.49),   0.30,   0.51,  False, "kinetic", True),
    ParamSpec("YEF_nom",     _g("YEF_nom",     0.49),   0.30,   0.51,  False, "kinetic", True),
    ParamSpec("Kn0_nom",     _g("Kn0_nom",     0.01),   1e-4,   1.0,   True,  "kinetic", True),
    ParamSpec("Kg0_nom",     _g("Kg0_nom",     7.5),    0.5,   50.0,   True,  "kinetic", True),
    ParamSpec("Kf0_nom",     _g("Kf0_nom",     7.5),    0.5,   50.0,   True,  "kinetic", True),
    ParamSpec("Kig0_nom",    _g("Kig0_nom",   55.0),    5.0,  500.0,   True,  "kinetic", True),
    ParamSpec("Kie0_nom",    _g("Kie0_nom",   40.0),    5.0,  500.0,   True,  "kinetic", True),
    ParamSpec("betaG0_nom",  _g("betaG0_nom", 0.225),   0.05,   1.0,   False, "kinetic", True),
    ParamSpec("betaF0_nom",  _g("betaF0_nom", 0.225),   0.05,   1.0,   False, "kinetic", True),
    ParamSpec("Kd0_nom",     _g("Kd0_nom",    4.4e-4),  1e-6,   1e-2,  True,  "kinetic", True),

    # === Grupo 2: Aroma soft constraints ===
    ParamSpec("PAIRWISE_PHI_STATIC",     _g("PAIRWISE_PHI_STATIC",     0.08),   1e-3, 0.5,  True,  "aroma", True),
    ParamSpec("PHI_ETHYL_ACETATE_STATIC",_g("PHI_ETHYL_ACETATE_STATIC",1e-4),   1e-6, 0.1,  True,  "aroma", True),   # soft coupling ethyl acetate (celda 42)

    # === Grupo 3: Turnover / N-reciclaje ===
    ParamSpec("TURNOVER_LAMBDA",     _g("TURNOVER_LAMBDA",    0.03),   1e-3,  0.2,   True,  "turnover", False),
    ParamSpec("K_NREC_UPTAKE",       _g("K_NREC_UPTAKE",      0.10),   1e-3,  1.0,   True,  "turnover", False),
    ParamSpec("Y_N_FROM_PROT",       _g("Y_N_FROM_PROT",      0.16),   0.10,  0.22,  False, "turnover", False),
    ParamSpec("K_AA_UPTAKE_GROWTH",  _g("K_AA_UPTAKE_GROWTH", 0.08),   1e-3,  1.0,   True,  "turnover", False),
]


# --- Herramientas de manipulacion ------------------------------------------------------
def active_specs() -> List[ParamSpec]:
    return [p for p in PARAM_SPECS if p.enabled]


def theta_to_physical(theta_unit: np.ndarray) -> Dict[str, float]:
    """
    Convierte vector en [0,1]^n (espacio de busqueda pymoo) a dict {name: valor_fisico}.
    Parametros log_scale se mapean: phys = 10**(log10(lb) + u * (log10(ub) - log10(lb))).
    """
    specs = active_specs()
    assert len(theta_unit) == len(specs), f"theta dim mismatch: {len(theta_unit)} vs {len(specs)}"
    out = {}
    for u, p in zip(theta_unit, specs):
        u = float(np.clip(u, 0.0, 1.0))
        if p.log_scale:
            val = 10 ** (np.log10(p.lb) + u * (np.log10(p.ub) - np.log10(p.lb)))
        else:
            val = p.lb + u * (p.ub - p.lb)
        out[p.name] = float(val)
    return out


def physical_to_theta(phys: Dict[str, float]) -> np.ndarray:
    specs = active_specs()
    u = np.zeros(len(specs))
    for i, p in enumerate(specs):
        v = float(phys.get(p.name, p.nominal))
        if p.log_scale:
            u[i] = (np.log10(v) - np.log10(p.lb)) / (np.log10(p.ub) - np.log10(p.lb))
        else:
            u[i] = (v - p.lb) / (p.ub - p.lb)
    return np.clip(u, 0.0, 1.0)


def apply_theta(theta_unit: np.ndarray, verbose: bool = False) -> Dict[str, float]:
    """
    Aplica theta (en espacio [0,1]) a las globales del notebook.
    Retorna dict con los valores fisicos aplicados (para logging).

    Si existe _sync_ester_soft_constraints() (registrada en el Modulo 3b), la invoca
    para actualizar las restricciones soft de esteres sobre `model` con los nuevos
    phi. Esto hace que las constraints se apliquen automaticamente en cada LP del
    forward.
    """
    phys = theta_to_physical(theta_unit)
    for name, val in phys.items():
        globals()[name] = float(val)
    if verbose:
        print("[apply_theta] Aplicado:")
        for k, v in phys.items():
            print(f"  {k:32s} = {v:.6g}")
    # Sync ester soft constraints if hook exists
    _hook = globals().get("_sync_ester_soft_constraints", None)
    if _hook is not None:
        try:
            _hook()
        except Exception as _e:
            if verbose:
                print(f"[apply_theta] WARNING: _sync_ester_soft_constraints fallo: {_e}")
    return phys


# --- Vector nominal (theta = u_nominal en [0,1]) --------------------------------------
theta_nominal_unit = physical_to_theta({p.name: p.nominal for p in active_specs()})

# --- Reporte resumen -------------------------------------------------------------------
print(f"[PARAM] Parametros activos: {len(active_specs())} / {len(PARAM_SPECS)}")
print(f"{'nombre':32s} {'grupo':10s} {'nominal':>12s} {'lb':>10s} {'ub':>10s} {'scale':>6s}")
for p in active_specs():
    print(f"{p.name:32s} {p.group:10s} {p.nominal:12.5g} {p.lb:10.3g} {p.ub:10.3g} "
          f"{'log' if p.log_scale else 'lin':>6s}")

N_PARAMS_ACTIVE = len(active_specs())
print(f"\nDimension del problema de optimizacion: n = {N_PARAMS_ACTIVE}")
