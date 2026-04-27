# Reporte: 02_dfba_hibrido_zenteno_henriques_clean

## Notebook generado

- Archivo: `02_dfba_hibrido_zenteno_henriques_clean.ipynb`.
- Objetivo: base limpia de dFBA COBRApy para fermentación vínica, sin PSO ni export MPCC.
- El archivo previo con el mismo nombre no era JSON válido y fue respaldado como `02_dfba_hibrido_zenteno_henriques_clean.ipynb.bak_invalid_20260421`.

## Bloques rescatados del notebook original

- Setup de rutas y carga de `yeast-GEM.xml`: celdas 3 y 5.
- IDs clave y ajuste de biomasa para anaerobiosis: celdas 7 y 9.
- Medio vínico anaeróbico y reglas de bounds: celda 11.
- Restricciones dinámicas tipo Zenteno: celda 17.
- Inclusión de aminoácidos aroma clave como fuentes de N cinéticas: leucina, valina, fenilalanina y metionina.
- Protein turnover, `N_rec` y GAM variable: celdas 38 y 50.
- Colocaciones ortogonales Radau IIA: celdas 51 y 52.
- Detector multifase simplificado: celda 53.

## Bloques descartados

- PSO y postprocesamiento de calibración: celdas 67-80.
- Diagnósticos extensos de aromas/redox/snapshots: celdas 40-49.
- Export MPCC, warm-start y metadata Julia/Notebook 2: celdas 55-64.
- Export para colaborador vía `exec(open(...).read())`: celdas 65-66.
- Variantes exploratorias duplicadas de SOA/FBA/pFBA: celdas 21-33, salvo lógica conceptual.

## Elementos tipo Zenteno

- Cinética de crecimiento dependiente de nitrógeno.
- Consumo de glucosa/fructosa con términos de crecimiento, fermentación y mantenimiento.
- Producción de etanol desde `betaG` y `betaF`.
- Inhibición por etanol y corrección térmica.
- Distribución de consumo de YAN sobre fuentes de N cinéticas.
- `N_gN_L` representa el YAN agregado no desagregado; `AA_leu_gN_L`, `AA_val_gN_L`, `AA_phe_gN_L` y `AA_met_gN_L` representan el bloque mínimo de aminoácidos aroma clave.

## Elementos tipo Henriques

- Función alternativa para hexosas con Michaelis-Menten e inhibición por etanol.
- Aromas fermentativos con ley directa `vP = kP * fermentation_hexose_Zenteno(T) * f_AA`.
- Caps suaves sobre exchanges de productos para mantener consistencia de orden de magnitud con el LP.
- Protein turnover activo en estacionaria/decay.
- GAM variable basado en composición macromolecular.
- Fases: lag, growth, n_limited_growth, stationary, decay.
- Corrección de fase para proteína: `lambda_turnover` y `vProt` solo operan en estacionaria/decay; `k_death_prot` solo opera en decay.
- Corrección de composición en decay: `Prot_g_L` y `Carb_g_L` pierden masa proporcional a la muerte de biomasa cuando `component_loss_tracks_biomass_death=True`, evitando inflación artificial de `Prot/X`.
- Separación de biomasa tipo Henriques: `XA_gDW_L` como biomasa activa/fermentativa, `XV_gDW_L` como biomasa viable tipo CFU y `XD_gDW_L` como biomasa muerta/inactiva; los flujos macroscópicos y aromas se escalan con `XA`.
- Bloque mínimo de aminoácidos aroma: leucina modula `isoamyl_acetate`; valina, fenilalanina y metionina quedan disponibles para isobutanol/isobutyl acetate, phenylethanol y methionol si se activan esos productos.
- Ensayo térmico aleatorio por punto de colocación: `temperature_profile_df` define escalones entre 15 y 25 °C que alimentan `dynamic_temperature`, las tasas Zenteno y los aromas directos Henriques.

## Funciones principales

- `validate_ids`
- `adjust_biomass_for_anaerobiosis`
- `apply_wine_medium`
- `quick_fba_check`
- `targeted_fva_check`
- `kinetic_limits_zenteno`
- `kinetic_limits_henriques_hexoses`
- `total_external_yan_gN_L`
- `amino_acid_uptake_caps`
- `aroma_aa_availability_factors`
- `aroma_product_factors`
- `biomass_pools`
- `viability_loss_rate`
- `product_rate_mass_action`
- `henriques_aroma_rates`
- `apply_henriques_product_caps`
- `compute_variable_gam`
- `apply_variable_gam`
- `compute_protein_turnover_bounds`
- `apply_protein_turnover`
- `infer_phase`
- `project_biomass_state`
- `solve_fba_at_state`
- `compute_rhs_from_fluxes`
- `run_dfba_collocation`

## Limitaciones actuales

- Solo se desagregó un bloque mínimo de aminoácidos aroma (`Leu`, `Val`, `Phe`, `Met`); el resto del YAN sigue agregado en `N_gN_L`.
- La producción de isoamyl acetate y ethyl acetate se integra como ODE directa tipo Henriques escalada por `XA_gDW_L`, pero `product_k` aún es nominal/no calibrado.
- La modulación por aminoácidos usa factores saturantes nominales (`aroma_aa_K_gN_L`, `aroma_aa_floor`); requiere calibración antes de interpretar magnitudes.
- La función de GAM variable replica la lógica global del notebook original, pero conviene validar coeficientes energéticos contra la reacción de biomasa específica.
- Si `Prot_g_L` o `Prot_frac_g_gDW` colapsan a cero tras esta corrección, eso indica que la fase estacionaria no está regenerando proteína vía `vProt` o que `turnover_lambda/k_nrec_uptake` están demasiado agresivos para los bounds disponibles.
- El smoke test es corto y no sustituye validación experimental.

## Próximos pasos sugeridos

1. Ejecutar el notebook completo y revisar diagnóstico FBA/FVA.
2. Ajustar `SMOKE_NFE`, `SMOKE_T_END_H` y luego correr una simulación nominal larga.
3. Probar `PARAMS["kinetic_form"] = "henriques_hexose"` y comparar contra Zenteno.
4. Calibrar `product_k` contra datos de isoamyl acetate y ethyl acetate antes de interpretar magnitudes.
5. Separar PSO/calibración en un notebook posterior.
