# Material para dFVB — Modelo dFBA Zenteno (yeast-GEM)

## Descripción General

Este directorio contiene todo el material necesario para implementar **dFVB (dynamic Flux Variability Analysis)** sobre el modelo de fermentación vínica basado en yeast-GEM.

El modelo corresponde a una **dFBA (dynamic Flux Balance Analysis)** de *Saccharomyces cerevisiae* con:
- Cinética tipo Zenteno (Monod + inhibición por etanol/glucosa)
- GAM variable (Growth-Associated Maintenance)
- Turnover de proteínas en fase estacionaria
- Seguimiento de aminoácidos individuales y producción de aromas (alcoholes superiores)

---

## Estructura de Archivos

```
Material_dFVB/
├── README.md                        ← Este archivo
├── export_model_for_dFVB.py         ← Script Python para generar los CSVs desde el notebook
│
├── import_model_data.m              ← MATLAB: script maestro de importación (cd a este dir)
├── load_params.m                    ← MATLAB: carga parámetros desde CSVs
├── ode_dfba_zenteno.m               ← MATLAB: función ODE del sistema dFBA
├── kinetic_limits_zenteno.m         ← MATLAB: cálculo de límites cinéticos (Zenteno)
├── compute_full_gam.m               ← MATLAB: cálculo de GAM variable
│
└── csv_exports/                     ← Generado al correr export_model_for_dFVB.py
    ├── S_matrix.csv                 ← Matriz estequiométrica densa (m x n)
    ├── S_matrix_sparse.csv          ← Matriz estequiométrica (formato triplete: row, col, value)
    ├── reaction_info.csv            ← ID, nombre, subsistema de cada reacción
    ├── metabolite_info.csv          ← ID, nombre, compartimento de cada metabolito
    ├── bounds_static.csv            ← Bounds estáticos completos (rxn_id, lb, ub, idx)
    ├── lb_vector.csv                ← Vector de cotas inferiores (n × 1, numérico puro)
    ├── ub_vector.csv                ← Vector de cotas superiores (n × 1, numérico puro)
    ├── bounds_dynamic_rules.csv     ← Reglas de bounds dinámicos por fase
    ├── objective_config.csv         ← Configuración de función objetivo por fase
    ├── phase_switching_criterion.csv← Criterio de cambio de fase
    ├── reaction_indices_of_interest.csv ← Índices de reacciones de interés (aromas, etc.)
    ├── initial_conditions.csv       ← Condiciones iniciales del estado
    ├── kinetic_parameters.csv       ← Parámetros cinéticos Zenteno
    ├── kinetic_param_set.txt        ← Nombre del set de parámetros usado
    ├── ode_state_definitions.csv    ← Definición del vector de estados ODE
    ├── gam_parameters.csv           ← Coeficientes GAM
    ├── nitrogen_config.csv          ← Mapa N_atoms por aminoácido
    ├── aroma_config.csv             ← Reacciones de intercambio de aromas
    ├── aa_precursor_config.csv      ← Reacciones de intercambio de aminoácidos precursores
    ├── turnover_parameters.csv      ← Parámetros de turnover de proteínas
    └── tracked_bounds_log.csv       ← Log de todas las modificaciones de bounds
```

---

## Paso a Paso para el Colaborador

### Paso 1: Generar los CSVs (Python)

> **Prerrequisito**: Ejecutar Notebook 1 (`01_preparacion_modelo_cobrapy.ipynb`) hasta la celda 46 inclusive, para que el modelo COBRApy y todas las variables estén en memoria.

Luego, ejecutar el contenido de `export_model_for_dFVB.py` como una celda nueva del notebook:

```python
exec(open("Material_dFVB/export_model_for_dFVB.py").read())
```

Esto genera la carpeta `Material_dFVB/csv_exports/` con todos los archivos listados arriba.

### Paso 2: Importar el modelo en MATLAB

```matlab
cd('ruta/a/Material_dFVB')
import_model_data    % script que carga todo al workspace
```

Esto deja en el workspace:
| Variable | Tipo | Descripción |
|----------|------|-------------|
| `S` | sparse matrix | Matriz estequiométrica [m × n] |
| `rxn_info` | table | Metadata de reacciones |
| `met_info` | table | Metadata de metabolitos |
| `lb`, `ub` | vector | Cotas de flujo |
| `rxn_indices` | table | Índices de reacciones categorizadas |
| `dyn_bounds` | table | Reglas de bounds dinámicos |
| `params` | struct | Todos los parámetros |
| `bounds_log` | table | Log de modificaciones de bounds |

### Paso 3: Usar las funciones MATLAB

```matlab
% Calcular límites cinéticos para un estado dado
lims = kinetic_limits_zenteno(G, F, E, N_free, X, params.kinetic);

% Actualizar GAM
fullGAM = compute_full_gam(P, R, C, Pbase, Rbase, Cbase, params.gam);

% ODE para integración
dydt = ode_dfba_zenteno(t, y, params, v_fba, phase);
```

---

## Descripción del Modelo dFBA

### Vector de Estados

El vector de estados `y` tiene la siguiente estructura:

| Índice | Estado | Unidad | Descripción |
|--------|--------|--------|-------------|
| 1 | X | gDW/L | Biomasa |
| 2 | N_free | gN/L | Nitrógeno libre (inorgánico) |
| 3 | G | g/L | Glucosa |
| 4 | F | g/L | Fructosa |
| 5 | E | g/L | Etanol |
| 6 | O2 | mmol/L | Oxígeno disuelto |
| 7 | Prot | g/L | Proteínas en biomasa |
| 8 | Carb | g/L | Carbohidratos en biomasa |
| 9..8+n_aa | AA_i | mmol/L | Aminoácidos individuales |
| 8+n_aa+1..end | Aroma_j | mmol/L | Alcoholes superiores |

### Fases de Operación

El modelo opera en **dos fases**:

1. **BIOMASS** (crecimiento): `N_total > 1e-3 gN/L`
   - Objetivo: maximizar crecimiento (`r_growth`)
   - Objetivo secundario: maximizar uptake de aminoácidos (ponderado por `w_aa`)
   - GAM variable según composición de biomasa

2. **TURNOVER_ATPM** (estacionaria/mantenimiento): `N_total ≤ 1e-3 gN/L`
   - Objetivo: maximizar ATP de mantenimiento (ATPM)
   - `r_growth = 0` (fijado con lb = ub = 0)
   - Proteínas se degradan a tasa `TURNOVER_LAMBDA`
   - N reciclado desde proteínas → pool de N inorgánico

### Criterio de Cambio de Fase

```
N_total = N_free + sum(N_atoms[aa_i] * MW_N * AA_i)    [gN/L]

Si N_total > 1e-3:  FASE = BIOMASS
Si N_total ≤ 1e-3:  FASE = TURNOVER_ATPM
```

### Cinética (Zenteno)

Parámetros predeterminados (set "paper2010"):

| Parámetro | Valor | Descripción |
|-----------|-------|-------------|
| MU0 | 0.18 h⁻¹ | Tasa máxima de crecimiento base |
| YXG | 1.6 gDW/gS | Rendimiento biomasa/sustrato glucosa |
| YXF | 1.6 gDW/gS | Rendimiento biomasa/sustrato fructosa |
| YXN | 19.69 gDW/gN | Rendimiento biomasa/nitrógeno |
| KG | 0.25 g/L | Constante Monod glucosa |
| KF | 0.25 g/L | Constante Monod fructosa |
| KN | 0.01 gN/L | Constante Monod nitrógeno |
| betaG0 | 0.225 L/g | Constante inhibición glucosa |
| betaF0 | 0.225 L/g | Constante inhibición fructosa |
| betaE0 | 0.008 L/g | Constante inhibición etanol |
| Ea_mu | 9729.0 K | Energía activación μ |
| Ea_Kg | -15000.0 K | Energía activación Kg |
| Ea_b | 7500.0 K | Energía activación β |
| T_ref | 303.15 K | Temperatura referencia |
| T_sim | 301.15 K | Temperatura simulación (28°C) |

La cinética incluye correcciones tipo Arrhenius para μ, Kg, y β:
```
mu_T  = MU0 * exp(Ea_mu * (1/T_ref - 1/T_sim))
Kg_T  = KG * exp(Ea_Kg * (1/T_ref - 1/T_sim))
b_T   = betaE0 * exp(Ea_b * (1/T_ref - 1/T_sim))
```

### GAM Variable

```
fullGAM = 30.49 + 16.965*(P/Pbase) + 1.638*(R/Rbase) + 5.210*Cfactor
Cfactor = max(0, (Cbase + Pbase - P - R) / Cbase)
```

Valores base: `Pbase ≈ 0.4612`, `Rbase ≈ 0.0600`, `Cbase ≈ 0.2873`

### Condiciones Iniciales (mosto tipo vino)

| Estado | Valor | Unidad |
|--------|-------|--------|
| X₀ | 0.5 | gDW/L |
| N₀ | 0.14 | gN/L |
| G₀ | 110.0 | g/L |
| F₀ | 110.0 | g/L |
| E₀ | 0.0 | g/L |
| O2₀ | 0.0 | mmol/L |
| Prot₀ | X₀ × 0.4612 | g/L |
| Carb₀ | X₀ × 0.2873 | g/L |
| AA₀ | Según concentración inicial en mosto | mmol/L |

### Modo Anaeróbico

El modelo opera en modo **anaeróbico**: intercambio de O₂ fijado con `lb = ub = 0`.

---

## Reacciones de Interés

### Aromas (Alcoholes Superiores)

Las reacciones de intercambio de aromas se encuentran en `aroma_config.csv`. Incluyen:
- Isobutanol (2-metil-1-propanol)
- Isoamyl alcohol (3-metil-1-butanol)
- 2-feniletanol
- Otros identificados por `found_aromas` en el notebook

### Precursores

Los aminoácidos precursores de aromas están identificados en `aa_config.csv`.

### Índices

`reaction_indices.csv` contiene los índices (0-based Python) de todas las reacciones categorizadas:
- `biomass`: reacción de crecimiento
- `atpm`: reacción de mantenimiento ATP
- `glucose_ex`: intercambio de glucosa
- `fructose_ex`: intercambio de fructosa
- `ethanol_ex`: intercambio de etanol
- `o2_ex`: intercambio de oxígeno
- `aroma_ex`: reacciones de intercambio de aromas
- `aa_ex`: reacciones de intercambio de aminoácidos

**Nota**: En MATLAB, sumar +1 a los índices (ya se hace en `import_model_data.m` con `rxn_indices.index_matlab`).

---

## Bounds Dinámicos

El archivo `dynamic_bounds.csv` documenta cómo cambian los bounds según la fase:

- **Fase BIOMASS**: Se aplican kinetic limits como upper bounds a las reacciones de uptake (glucosa, fructosa, N)
- **Fase TURNOVER_ATPM**: Se fija `r_growth = 0`, se desbloquea ATPM, se aplica reciclaje de N
- **Leaks cerrados**: Ciertas reacciones de fuga (`LEAK_CLOSE_IDS`) se fijan a 0 permanentemente
- **Sinks alternativos**: Ciertas reacciones sink tienen capacidad limitada (`ALT_SINK_CAPS`)

---

## Notas sobre el Modelo Base

- **Modelo**: yeast-GEM v8.x (archivo `yeast-GEM.xml`)
- **~7700 reacciones**, **~3800 metabolitos**
- Modificaciones aplicadas durante la preparación en Notebook 1 (Python/COBRApy):
  - Cierre de leaks
  - Limitación de sinks alternativos
  - Configuración anaeróbica
  - Bounds de médium (glucosa, fructosa, nitrógeno, aminoácidos)
  - Peso variable de ATPM en objetivo balanceado

Todas las modificaciones están documentadas en `bounds_log.csv`.

---

## Contacto

Para dudas sobre la implementación del modelo, referirse al notebook fuente:
`01_preparacion_modelo_cobrapy.ipynb`
