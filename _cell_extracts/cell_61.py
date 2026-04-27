# ======================================================================================
# MODULO 1 - Carga de datos experimentales (Excel multi-hoja)
# ======================================================================================
import numpy as np
import pandas as pd
from pathlib import Path

# --- Configuracion ---------------------------------------------------------------------
EXP_DATA_DIR  = Path("data_experimental")
EXP_XLSX_PATH = EXP_DATA_DIR / "Calibration_data_vl3.xlsx"

# Conversion Viability (10^6 cell/mL) -> gDW/L: 3.7e-12 g/cell * 1e6 cell/mL * 1e3 mL/L
VIABILITY_TO_GDW_L = 0.072

# Fallback sintetico si no hay Excel (delegado a celda 65 / Modulo 3c).
USE_SYNTHETIC_FROM_FORWARD = True
SYNTH_N_POINTS   = 15
SYNTH_NOISE_FRAC = 0.05
SYNTH_NOISE_FLOOR = 1e-3
SYNTH_SEED       = 42

# Sigma heuristica (fase 1): fraccion relativa + piso por observable.
SIGMA_FRAC = 0.05
SIGMA_FLOOR = {
    "X_gDW_L":              1e-3,
    "G_g_L":                2e-1,
    "F_g_L":                2e-1,
    "E_g_L":                2e-1,
    "N_gN_L":               5e-4,
    "isoamyl_acetate_mg_L": 1e-3,
    "ethyl_acetate_mg_L":   5e-2,
}

# Observables a ajustar (mismos nombres que la version sintetica).
TARGET_OBSERVABLES = [
    "X_gDW_L",
    "G_g_L",
    "F_g_L",
    "E_g_L",
    "N_gN_L",
    "isoamyl_acetate_mg_L",
    "ethyl_acetate_mg_L",
]

OBSERVABLE_WEIGHTS = {v: 1.0 for v in TARGET_OBSERVABLES}

# Mapa columna Excel -> observable + transformacion.
COLUMN_MAP = [
    ("Viability",              "X_gDW_L",              lambda v: v * VIABILITY_TO_GDW_L),
    ("GLUCOSE",                "G_g_L",                lambda v: v),
    ("FRUCTOSE",               "F_g_L",                lambda v: v),
    ("ETANOL",                 "E_g_L",                lambda v: v),
    ("YAN",                    "N_gN_L",               lambda v: v / 1000.0),
    ("isoamil_acetate_total",  "isoamyl_acetate_mg_L", lambda v: v),
    ("Ethyl_Acetate_total",    "ethyl_acetate_mg_L",   lambda v: v),
]


def _sigma_for(obs_name, value):
    floor = SIGMA_FLOOR.get(obs_name, 1e-6)
    return max(SIGMA_FRAC * abs(float(value)), floor)


def load_experimental_data_excel(path=EXP_XLSX_PATH):
    """
    Lee Excel multi-hoja. Cada hoja = ensayo independiente (dataset_id = nombre hoja).
    Retorna (exp_df long-form, experiments_meta dict).
    """
    path = Path(path)
    if not path.exists():
        print(f"[DATOS] {path} no encontrado.")
        return None, None

    xl = pd.ExcelFile(path)
    rows = []
    meta = {}

    for sheet in xl.sheet_names:
        df = xl.parse(sheet).sort_values("t").reset_index(drop=True)
        ds_id = str(sheet)

        # --- observaciones en long-form ---
        for col_src, obs_name, transform in COLUMN_MAP:
            if col_src not in df.columns:
                continue
            for _, r in df.iterrows():
                v = r[col_src]
                if pd.isna(v):
                    continue
                y = float(transform(v))
                rows.append({
                    "dataset_id": ds_id,
                    "t_h":        float(r["t"]),
                    "variable":   obs_name,
                    "value":      y,
                    "sigma":      float(_sigma_for(obs_name, y)),
                })

        # --- condiciones iniciales (fila t=0) ---
        r0 = df.iloc[0]

        def _get(col, default, tr=lambda x: x):
            v = r0.get(col, np.nan)
            return float(tr(v)) if pd.notna(v) else default

        x0 = _get("Viability", default=None, tr=lambda v: v * VIABILITY_TO_GDW_L)
        if x0 is None:
            # fallback a Peso Seco si Viability t=0 es NaN
            x0 = _get("Peso Seco", default=0.08)

        initial_conditions = {
            "X_gDW_L": x0,
            "G_g_L":   _get("GLUCOSE",  default=110.0),
            "F_g_L":   _get("FRUCTOSE", default=110.0),
            "E_g_L":   _get("ETANOL",   default=0.0),
            "N_gN_L":  _get("YAN",      default=0.14, tr=lambda v: v / 1000.0),
        }

        # Aromas IC (mg/L) desde fila t=0.
        aroma_ic_mg_L = {
            "isoamyl_acetate_mg_L": _get("isoamil_acetate_total", default=0.0),
            "ethyl_acetate_mg_L":   _get("Ethyl_Acetate_total",   default=0.0),
        }

        # --- perfil de temperatura (C -> K) ---
        temp_df = (df[["t", "temperatura"]].dropna()
                     .rename(columns={"t": "t_h", "temperatura": "T_C"}))
        temp_df["T_K"] = temp_df["T_C"] + 273.15
        T_profile = temp_df[["t_h", "T_K"]].reset_index(drop=True)

        # --- pulsos de N ---
        pulses = df[df["pulso_nut"].notna()][["t", "pulso_nut"]]
        pulse_schedule = [
            {"t_h": float(r["t"]), "amount_gN_L": float(r["pulso_nut"]) / 1000.0}
            for _, r in pulses.iterrows()
        ]

        meta[ds_id] = {
            "dataset_id":          ds_id,
            "t_end":               float(df["t"].max()),
            "n_points":            int(len(df)),
            "initial_conditions":  initial_conditions,
            "aroma_ic_mg_L":       aroma_ic_mg_L,
            "T_profile":           T_profile,
            "pulse_schedule":      pulse_schedule,
        }

    exp_df = pd.DataFrame(rows)
    exp_df = exp_df[exp_df["variable"].isin(TARGET_OBSERVABLES)].reset_index(drop=True)
    return exp_df, meta


# --- Carga efectiva --------------------------------------------------------------------
EXP_DATA_DIR.mkdir(exist_ok=True)
exp_df, experiments_meta = load_experimental_data_excel()
if experiments_meta is None:
    experiments_meta = {}

# ──────────────────────────────────────────────────────────────────────────────
# Seleccion manual de hojas para calibracion (editable por el usuario).
# Vacio = todas las hojas detectadas; de lo contrario, lista de dataset_id a usar.
# ──────────────────────────────────────────────────────────────────────────────
DATASET_WHITELIST = []   # ej: ["25026", "25085"]

_all_ids = list(experiments_meta.keys())
if DATASET_WHITELIST:
    _unknown = [d for d in DATASET_WHITELIST if d not in experiments_meta]
    if _unknown:
        raise ValueError(f"DATASET_WHITELIST contiene ids no presentes: {_unknown}. "
                         f"Disponibles: {_all_ids}")
    SELECTED_DATASETS = list(DATASET_WHITELIST)
else:
    SELECTED_DATASETS = list(_all_ids)
DATASET_PRIMARY = SELECTED_DATASETS[0] if SELECTED_DATASETS else None


# --- Helpers entradas de proceso por experimento (Fases 2, 3, 4) ----------------------
_PROC_GLOBAL_KEYS = (
    # IC (Fase 2)
    "X0_init", "G0_init", "F0_init", "E0_init", "N0",
    # T profile (Fase 3) y pulsos (Fase 4)
    "T_PROFILE_DF", "PULSE_SCHEDULE",
    # Aromas IC (mg/L) desde t=0
    "AROMA_IC_MG_L",
)


def set_process_inputs(ds_id):
    """
    Inyecta IC + T-profile + pulsos del ensayo `ds_id` como globales que lee el solver
    dFBA (celda 44). Retorna snapshot previo para restaurar.
    """
    if ds_id not in experiments_meta:
        raise KeyError(f"dataset {ds_id!r} no esta en experiments_meta")
    meta = experiments_meta[ds_id]
    ic   = meta["initial_conditions"]

    snap = {k: globals().get(k, None) for k in _PROC_GLOBAL_KEYS}

    globals()["X0_init"]       = float(ic["X_gDW_L"])
    globals()["G0_init"]       = float(ic["G_g_L"])
    globals()["F0_init"]       = float(ic["F_g_L"])
    globals()["E0_init"]       = float(ic["E_g_L"])
    globals()["N0"]            = float(ic["N_gN_L"])
    globals()["T_PROFILE_DF"]  = meta["T_profile"]
    globals()["PULSE_SCHEDULE"] = list(meta["pulse_schedule"])
    globals()["AROMA_IC_MG_L"] = dict(meta.get("aroma_ic_mg_L", {}))
    return snap


def restore_process_inputs(snapshot):
    for k, v in snapshot.items():
        if v is None:
            globals().pop(k, None)
        else:
            globals()[k] = v


# Aliases retro-compatibles con nombres anteriores
set_ic_from_experiment = set_process_inputs
restore_ic_globals     = restore_process_inputs


# --- Reporte ---------------------------------------------------------------------------
if exp_df is not None and len(exp_df):
    print(f"[DATOS] {len(exp_df)} observaciones sobre {exp_df['dataset_id'].nunique()} "
          f"ensayos, {exp_df['variable'].nunique()} variables objetivo.")
    print("\n[DATOS] n puntos por ensayo x variable:")
    print(exp_df.groupby(["dataset_id", "variable"]).size().unstack(fill_value=0))

    print("\n[DATOS] Metadata por ensayo:")
    for ds, m in experiments_meta.items():
        pulses = len(m["pulse_schedule"])
        Tmin = m["T_profile"]["T_K"].min() - 273.15 if len(m["T_profile"]) else float("nan")
        Tmax = m["T_profile"]["T_K"].max() - 273.15 if len(m["T_profile"]) else float("nan")
        ic = m["initial_conditions"]
        print(f"  {ds}: t_end={m['t_end']:6.0f} h | T={Tmin:5.1f}-{Tmax:5.1f} C | "
              f"pulsos={pulses} | IC X={ic['X_gDW_L']:.3f} N={ic['N_gN_L']:.4f} "
              f"G={ic['G_g_L']:.1f} F={ic['F_g_L']:.1f} E={ic['E_g_L']:.2f}")
    print(f"\n[DATOS] SELECTED_DATASETS = {SELECTED_DATASETS}")
    print(f"[DATOS] DATASET_PRIMARY    = {DATASET_PRIMARY}")
else:
    print("[DATOS] exp_df no cargado (Excel ausente). Se intentara fallback sintetico "
          "en celda 65 si USE_SYNTHETIC_FROM_FORWARD=True.")
    exp_df = None
