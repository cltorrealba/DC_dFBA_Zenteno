# Variante balanceada basada en celda 45:
# 1) Mantener crecimiento maximo (mu*) en fase BIOMASS.
# 2) ATPM liberado en fase BIOMASS.
# 3) Objetivo mixto desactivado por defecto para correr equivalente a celda 45.

# -----------------------------
# Parametros de la variante
# -----------------------------
BAL_GROWTH_LOCK_FRAC = 1.0   # 1.0 = exige mu* (usa 0.999 si quieres robustez numerica)
BAL_W_ATPM_GROWTH    = 0.0   # 0.0 = sin sesgo ATPM (modo equivalente a celda 45)
BAL_W_AA_GROWTH      = 0.0   # 0.0 = sin sesgo uptake AA
BAL_PFBA_MIX         = True  # si activas mix, aplica pFBA al objetivo mixto

# ATPM libre en crecimiento (en vez de fijo lb=ub=0.7)
BAL_ATPM_GROWTH_LB   = 0.7
BAL_ATPM_GROWTH_UB   = 1000.0


def _fba_at_state_colloc_balanced(
    yk,
    aerobic_mode: bool,
    growth_lock_frac: float = BAL_GROWTH_LOCK_FRAC,
    w_atpm_growth: float = BAL_W_ATPM_GROWTH,
    w_aa_growth: float = BAL_W_AA_GROWTH,
    use_pfba_mix: bool = BAL_PFBA_MIX,
):
    """
    Igual que _fba_at_state_colloc, pero en BIOMASS aplica:
      (A) max mu -> obtiene mu*
      (B) fija mu >= growth_lock_frac * mu*
      (C) opcional: max objetivo mixto (ATPM + AA uptake) dentro del mismo crecimiento

    Nota: con w_atpm_growth=0 y w_aa_growth=0, se mantiene modo equivalente a celda 45
    (sin sesgo secundario), pero con ATPM liberado en fase de crecimiento.
    """
    yk = np.maximum(yk, 0.0)
    cX       = yk[IDX_X]
    cN_free  = yk[IDX_NFREE]
    cG       = yk[IDX_G]
    cF       = yk[IDX_F]
    cE       = yk[IDX_E]
    cO2      = yk[IDX_O2]
    cProt    = yk[IDX_PROT]
    cCarb    = yk[IDX_CARB]
    cNrec    = yk[IDX_NREC]
    aa_state = {ak: float(yk[IDX_AA0 + j]) for j, ak in enumerate(aa_keys)}

    aa_n_sum = MW_N * sum(aa_state.values())
    cN_total = cN_free + aa_n_sum

    P_frac   = cProt / max(cX, EPS_C)
    C_frac   = cCarb / max(cX, EPS_C)
    full_gam = _compute_full_gam(P_frac, RNA_FRAC, C_frac,
                                 Pbase_global, Rbase_global, Cbase_global)

    lim_glu, lim_fru, lim_eth, lim_obj, lim_n_map = _kinetic_limits_fn_c(
        cX, cN_free, cG, cF, cE, cO2)

    with model as mtmp:
        rb = {r.id: r for r in mtmp.reactions}

        for rid in _tracked_ids:
            if rid in rb:
                rb[rid].lower_bound = _base_lb_nrec[rid]
                rb[rid].upper_bound = _base_ub_nrec[rid]

        if GLU_ID in rb:
            rb[GLU_ID].lower_bound = max(rb[GLU_ID].lower_bound, -lim_glu)
            rb[GLU_ID].upper_bound = min(rb[GLU_ID].upper_bound, 0.0)
        if FRU_ID in rb:
            rb[FRU_ID].lower_bound = max(rb[FRU_ID].lower_bound, -lim_fru)
            rb[FRU_ID].upper_bound = min(rb[FRU_ID].upper_bound, 0.0)

        for rid in KINETIC_N_SOURCE_IDS:
            if rid in rb:
                Li = float(lim_n_map.get(rid, 0.0))
                rb[rid].lower_bound = max(rb[rid].lower_bound, -Li)
                rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)

        if aerobic_mode:
            lim_o2 = O2_VMAX_UPTAKE * (cO2 / (cO2 + O2_KO_G_L + EPS_C))
            lim_o2 = 0.0 if cO2 <= O2_DEPLETION_THRESHOLD else max(0.0, lim_o2)
            rb[O2_ID].lower_bound = -lim_o2
            rb[O2_ID].upper_bound = 0.0
        else:
            rb[O2_ID].lower_bound = 0.0
            rb[O2_ID].upper_bound = 0.0

        if APPLY_PRODUCT_CAPS:
            rb[ETH_ID].upper_bound = min(rb[ETH_ID].upper_bound, lim_eth)
            rb[OBJ_ID].upper_bound = min(rb[OBJ_ID].upper_bound, lim_obj)
        else:
            rb[ETH_ID].upper_bound = max(rb[ETH_ID].upper_bound, 1000.0)
            rb[OBJ_ID].upper_bound = max(rb[OBJ_ID].upper_bound, 1000.0)

        _apply_variable_gam(mtmp, full_gam)

        if cN_total > N_TOTAL_DEPLETION_THRESHOLD:
            # Fase crecimiento
            for ak, rid in found_precursors.items():
                if rid in rb:
                    q_i = max(0.0, K_AA_UPTAKE_GROWTH * aa_state[ak] / max(cX, EPS_C))
                    rb[rid].lower_bound = -q_i
                    rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)

            # ATPM libre en crecimiento
            if ATPM_ID in rb:
                rb[ATPM_ID].lower_bound = max(rb[ATPM_ID].lower_bound, BAL_ATPM_GROWTH_LB)
                rb[ATPM_ID].upper_bound = max(rb[ATPM_ID].upper_bound, BAL_ATPM_GROWTH_UB)

            # Paso A: max crecimiento
            mtmp.objective = OBJ_ID
            mtmp.objective_direction = "max"
            try:
                sol_mu = mtmp.optimize()
            except Exception:
                return None, "FAIL"
            if getattr(sol_mu, "status", "") != "optimal":
                return None, "NO_OPT"

            mu_star = max(0.0, float(sol_mu.fluxes.get(OBJ_ID, 0.0)))

            # Paso B: bloquear crecimiento maximo
            mu_lb = max(0.0, float(growth_lock_frac) * mu_star)
            rb[OBJ_ID].lower_bound = max(rb[OBJ_ID].lower_bound, mu_lb)

            # Paso C (opcional): objetivo mixto ATPM + uptake de AA
            mix_obj = {}
            if ATPM_ID in rb and w_atpm_growth > 0.0:
                atpm_ref = max(abs(float(sol_mu.fluxes.get(ATPM_ID, 0.0))), 1e-9)
                mix_obj[mtmp.reactions.get_by_id(ATPM_ID)] = float(w_atpm_growth) / atpm_ref

            if w_aa_growth > 0.0:
                uptake_ref = max(1.0, float(len(found_precursors)))
                for rid in found_precursors.values():
                    if rid in rb:
                        # uptake es flujo negativo en EX_*; coef negativo favorece mayor uptake
                        mix_obj[mtmp.reactions.get_by_id(rid)] = -float(w_aa_growth) / uptake_ref

            if mix_obj:
                mtmp.objective = mix_obj
                mtmp.objective_direction = "max"
                try:
                    sol = cobra.flux_analysis.pfba(mtmp, fraction_of_optimum=1.0) if use_pfba_mix else mtmp.optimize()
                except Exception:
                    return None, "FAIL"
            else:
                # Modo equivalente a celda 45: conservar solucion de max crecimiento
                sol = sol_mu

            mode = "BIOMASS"

        else:
            # Fase estacionaria (igual a celda 45)
            for rid in KINETIC_N_SOURCE_IDS:
                if rid in rb:
                    rb[rid].lower_bound = max(rb[rid].lower_bound, 0.0)
                    rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)
            for rid in N_VITAMIN_EX_IDS_C:
                if rid in rb:
                    rb[rid].lower_bound = max(rb[rid].lower_bound, 0.0)
                    rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)

            rb[OBJ_ID].lower_bound = 0.0
            rb[OBJ_ID].upper_bound = 0.0

            q_Nrec = max(0.0, K_NREC_UPTAKE * cNrec / max(cX, EPS_C))
            q_AA   = q_Nrec / MW_N
            for ak, rid in found_precursors.items():
                if rid in rb:
                    rb[rid].lower_bound = -max(0.0, q_AA * AA_ALPHA[ak])
                    rb[rid].upper_bound = min(rb[rid].upper_bound, 0.0)

            vprot_star = 0.0
            try:
                mtmp.objective = PROT_RXN_ID
                mtmp.objective_direction = "max"
                sol_p = mtmp.optimize()
                if getattr(sol_p, "status", "") == "optimal":
                    vprot_star = max(0.0, float(sol_p.fluxes.get(PROT_RXN_ID, 0.0)))
            except Exception:
                vprot_star = 0.0

            if vprot_star > 1e-10:
                rb[PROT_RXN_ID].lower_bound = max(rb[PROT_RXN_ID].lower_bound,
                                                  VPROT_FLOOR_FRAC * vprot_star)

            rb[ATPM_ID].lower_bound = max(rb[ATPM_ID].lower_bound, ATPM_LB_NO_GROWTH)
            rb[ATPM_ID].upper_bound = ATPM_UB_NO_GROWTH

            mtmp.objective = ATPM_ID
            mtmp.objective_direction = "max"
            try:
                sol = cobra.flux_analysis.pfba(mtmp, fraction_of_optimum=1.0)
            except Exception:
                return None, "FAIL"
            mode = "TURNOVER_ATPM"

        if getattr(sol, "status", "optimal") != "optimal":
            return None, "NO_OPT"

        fluxes = dict(
            v_obj    = float(sol.fluxes.get(OBJ_ID,      0.0)),
            v_glu    = float(sol.fluxes.get(GLU_ID,      0.0)),
            v_fru    = float(sol.fluxes.get(FRU_ID,      0.0)),
            v_eth    = float(sol.fluxes.get(ETH_ID,      0.0)),
            v_o2     = float(sol.fluxes.get(O2_ID,       0.0)),
            v_atpm   = float(sol.fluxes.get(ATPM_ID,     0.0)),
            v_prot   = float(sol.fluxes.get(PROT_RXN_ID, 0.0)),
            full_gam = full_gam,
        )
        vn_eff = 0.0
        if mode == "BIOMASS":
            for rid in KINETIC_N_SOURCE_IDS:
                vi = float(sol.fluxes.get(rid, 0.0))
                vn_eff += max(0.0, -vi) * N_atoms_map.get(rid, 1.0) * MW_N
        fluxes["vn_eff"] = vn_eff

        for ak in aa_keys:
            fluxes[f"v_aa_{ak}"] = float(sol.fluxes.get(found_precursors[ak], 0.0))
        for ak_a in aroma_keys:
            fluxes[f"v_aroma_{ak_a}"] = float(sol.fluxes.get(found_aromas[ak_a], 0.0))

        # Fluxes extra por rid (acetate esters que no estan en found_aromas).
        for _rid in globals().get("EXTRA_AROMA_RIDS", []):
            fluxes[f"v_aroma_{_rid}"] = float(sol.fluxes.get(_rid, 0.0))

    return fluxes, mode


def _run_dfba_colloc_balanced(
    growth_lock_frac=BAL_GROWTH_LOCK_FRAC,
    w_atpm_growth=BAL_W_ATPM_GROWTH,
    w_aa_growth=BAL_W_AA_GROWTH,
    use_pfba_mix=BAL_PFBA_MIX,
    **kwargs,
):
    """
    Wrapper de _run_dfba_colloc que reutiliza toda la logica de celda 45,
    cambiando solo el solver local de cada punto por la variante balanceada.
    """
    _original = _fba_at_state_colloc

    def _balanced_adapter(yk, aerobic_mode):
        return _fba_at_state_colloc_balanced(
            yk,
            aerobic_mode,
            growth_lock_frac=growth_lock_frac,
            w_atpm_growth=w_atpm_growth,
            w_aa_growth=w_aa_growth,
            use_pfba_mix=use_pfba_mix,
        )

    globals()["_fba_at_state_colloc"] = _balanced_adapter
    try:
        return _run_dfba_colloc(**kwargs)
    finally:
        globals()["_fba_at_state_colloc"] = _original


print("Implementacion balanceada lista: _run_dfba_colloc_balanced")
print("\n" + "=" * 80)
print("Escenario: Anaerobico + N_rec + Radau IIA + ATPM_liberado (sin mix)  (nfe=18, ncp=3, n_iter=2)")
print("=" * 80)

res_colloc_bal = _run_dfba_colloc_balanced(
    scenario_name="Anaerobico + N_rec + Radau IIA + ATPM_liberado (sin mix)",
    aerobic_mode=False,
    o2_init_g_l=0.0,
    nfe=18,
    ncp=3,
    t_end=72.0,
    n_iter=2,
    growth_lock_frac=1.0,
    w_atpm_growth=0.0,
    w_aa_growth=0.0,
    use_pfba_mix=True,
)

# Exportar warm-start balanceado en formato NB2
_C_ws_bal = res_colloc_bal["C"]
_CDOT_ws_bal = res_colloc_bal["CDOT"]
_hm_ws_bal = res_colloc_bal["hm"]
_nfe_ws_bal = res_colloc_bal["nfe"]
_ncp_ws_bal = res_colloc_bal["ncp"]

_rows_c_bal = []
_rows_cdot_bal = []
for _s in range(n_state):
    for _i in range(_nfe_ws_bal):
        for _j in range(_ncp_ws_bal):
            _rows_c_bal.append({"state_idx": _s, "fe": _i, "cp": _j, "value": float(_C_ws_bal[_s, _i, _j])})
            _rows_cdot_bal.append({"state_idx": _s, "fe": _i, "cp": _j, "value": float(_CDOT_ws_bal[_s, _i, _j])})

_ws_c_file_bal = OUT_DIR / "warm_start_states_robust_colloc_balanced.csv"
_ws_cdot_file_bal = OUT_DIR / "warm_start_cdot_robust_colloc_balanced.csv"
_ws_h_file_bal = OUT_DIR / "warm_start_h_robust_fe_balanced.csv"

pd.DataFrame(_rows_c_bal).to_csv(_ws_c_file_bal, index=False)
pd.DataFrame(_rows_cdot_bal).to_csv(_ws_cdot_file_bal, index=False)
np.savetxt(_ws_h_file_bal, _hm_ws_bal, delimiter=",")

# Actualizar referencias globales para celdas downstream
states_df_phase = res_colloc_bal["states_df"].copy()
flux_df_phase = res_colloc_bal["flux_df"].copy()

print("\nWarm-start balanceado exportado (listo para NB2):")
print(f"  {_ws_c_file_bal}")
print(f"  {_ws_cdot_file_bal}")
print(f"  {_ws_h_file_bal}")
print(f"  ODE gap = {res_colloc_bal['max_ode_gap']:.2e}")
print("[INFO] states_df_phase y flux_df_phase actualizados con trayectoria balanceada.")