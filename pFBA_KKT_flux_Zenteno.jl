# =============================================================================
# pFBA_KKT_flux_Zenteno_O2minimal.jl
# -----------------------------------------------------------------------------
# - O2 como estado dinámico (balance), sin restricción cinética explícita.
# - Agrega flux constraint (paper-inspired) ENTRE elementos:
#       |v(fe) - v(fe-1)| <= eps_flux
# =============================================================================

function pFBA_KKT_flux_Zenteno_O2minimal(
    c0;
    eps_flux::Float64 = 0.0,
    apply_product_caps::Bool = false
)

  # Estados cineticos de Zenteno: X,N,G,F,E,O2
  @assert nc == 6 "Este script asume nc=6 (X,N,G,F,E,O2)."
  @assert length(c0) == 6 "c0 debe tener 6 estados: [X0,N0,G0,F0,E0,O20]."
  o2_tol = try parse(Float64, get(ENV, "O2_TOL", "1e-8")) catch; 1e-8 end

  # --- Radau IIA 3 points
  colmat = [0.19681547722366  -0.06553542585020  0.02377097434822;
            0.39442431473909   0.29207341166523 -0.04154875212600;
            0.37640306270047   0.51248582618842  0.11111111111111]
  radau  = [0.15505, 0.64495, 1.00000]

  # --- Model
  m = Model(optimizer_with_attributes(
      Ipopt.Optimizer,
      "linear_solver" => "mumps",
      "print_level" => 5,
      "tol" => 1e-4,
      "acceptable_tol" => 1e-2,
      "acceptable_iter" => 12,
      "acceptable_constr_viol_tol" => 1e-3,
      "acceptable_compl_inf_tol" => 1e-3,
      "max_iter" => 100,
      "constr_viol_tol" => 1e-5,
      "compl_inf_tol" => 1e-4,
      "mu_strategy" => "adaptive",
      "mu_init" => 1e-2,
      "mu_min" => 1e-8,
    #   "nlp_scaling_method" => "equilibration-based",
      "nlp_scaling_max_gradient" => 100.0,
      "obj_scaling_factor" => 1e-4,
      "bound_relax_factor" => 1e-8,
      "warm_start_init_point" => "yes"
  ))

# m = Model(optimizer_with_attributes(Ipopt.Optimizer, "warm_start_init_point" => "yes", "print_level" => 5, "linear_solver" => "mumps"))


  # --- Variables
  @variables(m, begin
      c[1:nc, 1:nfe, 1:ncp] >= 0.0
      cdot[1:nc, 1:nfe, 1:ncp]
      v[1:nv, 1:nfe]
      lambda_[1:nm, 1:nfe]
      alpha_U[1:nv, 1:nfe] >= 0.0
      alpha_L[1:nv, 1:nfe] <= 0.0
      alpha_upt[1:n_up, 1:nfe] <= 0.0
      FO_U[1:nv, 1:nfe]
      FO_L[1:nv, 1:nfe]
      FO_upt[1:n_up, 1:nfe]
      alpha_prod[1:n_prod, 1:nfe] >= 0.0
      FO_prod[1:n_prod, 1:nfe]
      hv[1:nfe] >= 0.0
  end)

  # --- Start values
  for i in 1:nfe
      set_start_value(hv[i], hm[i])
      set_start_value(v[o2,i], 0.0)
      for j in 1:ncp
          set_start_value(c[6,i,j], c0[6] / cs[6])
      end
  end

  # --- Scale initial conditions
  c0s = copy(c0)
  for i in 1:nc
      c0s[i] = c0s[i] / cs[i]
  end

  # --- Register custom functions
  JuMP.register(m, :dynamic_temperature, 1, dynamic_temperature; autodiff=true)
  JuMP.register(m, :death_rate_T, 2, death_rate_T; autodiff=true)
  JuMP.register(m, :smooth_injection, 4, smooth_injection; autodiff=true)

  # --- Time mapping (hv variable)
  @NLexpression(m, t_fe[i=1:nfe], sum(hv[k] for k in 1:i-1))
  @NLexpression(m, t_loc[i=1:nfe, j=1:ncp], t_fe[i] + radau[j]*hv[i])

  # --- Temperature / Arrhenius
  @NLexpression(m, T_loc[i=1:nfe, j=1:ncp], dynamic_temperature(t_loc[i,j]))
  @NLexpressions(m, begin
      # c variables are scaled; define physical concentrations explicitly.
      cX[i=1:nfe, j=1:ncp], c[1,i,j] * cs[1]
      cN[i=1:nfe, j=1:ncp], c[2,i,j] * cs[2]
      cG[i=1:nfe, j=1:ncp], c[3,i,j] * cs[3]
      cF[i=1:nfe, j=1:ncp], c[4,i,j] * cs[4]
      cE[i=1:nfe, j=1:ncp], c[5,i,j] * cs[5]
  end)
  @NLexpression(m, mu_T_ij[i=1:nfe, j=1:ncp],
      exp(59453.0 * (T_loc[i,j] - 300.0) / (300.0 * R * T_loc[i,j])))
  @NLexpression(m, Kg_T_ij[i=1:nfe, j=1:ncp],
      exp(46055.0 * (T_loc[i,j] - 293.15) / (293.15 * R * T_loc[i,j])))
  @NLexpression(m, b_T_ij[i=1:nfe, j=1:ncp],
      exp(11000.0 * (T_loc[i,j] - 296.15) / (296.15 * R * T_loc[i,j])))
  @NLexpression(m, mrate_ij[i=1:nfe, j=1:ncp],
      MRATE_0 * exp(37681.0 * (T_loc[i,j] - 293.30) / (293.30 * R * T_loc[i,j])))

  # --- Zenteno kinetics
  @NLexpression(m, mu_j[i=1:nfe, j=1:ncp],
      MU0 * mu_T_ij[i,j] * (cN[i,j] / (cN[i,j] + Kn0_nom * Kg_T_ij[i,j] + EPS)))

  @NLexpression(m, betaG_j[i=1:nfe, j=1:ncp],
      betaG0_nom * b_T_ij[i,j] *
      (cG[i,j] / (cG[i,j] + Kg0_nom * Kg_T_ij[i,j] + EPS)) *
      (Kie0_nom * Kg_T_ij[i,j] / (cE[i,j] + Kie0_nom * Kg_T_ij[i,j] + EPS)))

  @NLexpression(m, betaF_j[i=1:nfe, j=1:ncp],
      betaF0_nom * b_T_ij[i,j] *
      (cF[i,j] / (cF[i,j] + Kf0_nom * Kg_T_ij[i,j] + EPS)) *
      (Kig0_nom * Kg_T_ij[i,j] / (cG[i,j] + Kig0_nom * Kg_T_ij[i,j] + EPS)) *
      (Kie0_nom * Kg_T_ij[i,j] / (cE[i,j] + Kie0_nom * Kg_T_ij[i,j] + EPS)))

  @NLexpression(m, Kd_j[i=1:nfe, j=1:ncp],
      death_rate_T(cE[i,j], T_loc[i,j]))

  # --- Nutrient injection
  @NLexpression(m, injection_rate[i=1:nfe, j=1:ncp],
      smooth_injection(t_loc[i,j], T_INJ_1, DOSE_1, WIDTH_1) +
      smooth_injection(t_loc[i,j], T_INJ_2, DOSE_2, WIDTH_2))

  # --- Dynamic macro limits
    PROD_CAP_ON = apply_product_caps ? 1.0 : 0.0
    PROD_CAP_BIG = 1e6
  @NLexpressions(m, begin
      vx[i=1:nfe, j=1:ncp],  mu_j[i,j]

      vg[i=1:nfe, j=1:ncp],  (mu_j[i,j] / YXG_nom +
                     betaG_j[i,j] / YEG +
                     mrate_ij[i,j] * (cG[i,j] / (cG[i,j] + cF[i,j] + EPS)))

      vf[i=1:nfe, j=1:ncp],  (mu_j[i,j] / YXF_nom +
                     betaF_j[i,j] / YEF +
                     mrate_ij[i,j] * (cF[i,j] / (cG[i,j] + cF[i,j] + EPS)))

      vn[i=1:nfe, j=1:ncp],  (mu_j[i,j] / YXN)

      L_uptake[k=1:n_up, i=1:nfe, j=1:ncp],
          IS_GLU[k] * (vg[i,j] / MW_GLU) +
          IS_FRU[k] * (vf[i,j] / MW_FRU) +
          IS_NIT[k] * vn[i,j] * N_frac[k]

      ve[i=1:nfe, j=1:ncp],  (betaG_j[i,j] + betaF_j[i,j]) / MW_ETH

      L_product[k=1:n_prod, i=1:nfe, j=1:ncp],
          IS_ETH_prod[k] * ve[i,j] + IS_OBJ_prod[k] * vx[i,j]

      # Si apply_product_caps=false, se relaja con un limite grande (sin capeo efectivo)
      L_product_cap[k=1:n_prod, i=1:nfe],
          L_product[k,i,ncp] + (1.0 - PROD_CAP_ON) * PROD_CAP_BIG
  end)

  # --- Objective
  @NLobjective(m, Min,
      sum(sum(-phi1*FO_L[mc,i] - phi3*FO_U[mc,i] for mc in 1:nv) +
          sum(phi2*FO_upt[k,i] for k in 1:n_up) +
          sum(-phi4*FO_prod[k,i] for k in 1:n_prod)
      for i in 1:nfe))

  # --- Linear constraints
  @constraints(m, begin
      # Collocation
      coll_c_0[l=1:nc, j=1:ncp],
          c[l,1,j] == c0s[l] + hv[1] * sum(colmat[j,k] * cdot[l,1,k] for k in 1:ncp)

      coll_c_n[l=1:nc, i=2:nfe, j=1:ncp],
          c[l,i,j] == c[l,i-1,ncp] + hv[i] * sum(colmat[j,k] * cdot[l,i,k] for k in 1:ncp)

      # FBA
      Sc[mc=1:nm, i=1:nfe],
          sum(S[mc,k] * v[k,i] * vs[k] for k in 1:nv) == 0

      v_UB[mc=1:nv, i=1:nfe],  v[mc,i]*vs[mc] - vub[mc] <= 0
      v_LB[mc=1:nv, i=1:nfe], -v[mc,i]*vs[mc] + vlb[mc] <= 0
      anaerobic_o2_pos[i=1:nfe],  v[o2,i]*vs[o2] <= o2_tol
      anaerobic_o2_neg[i=1:nfe], -v[o2,i]*vs[o2] <= o2_tol

      # Mesh
      MFE1, sum(hv[i] for i in 1:nfe) == th
      MFE3[i=1:nfe], hv[i] >= 0.0
      MFE4[i=1:nfe], hv[i] >= (1.0 - var_h) * hm[1]
      MFE5[i=1:nfe], hv[i] <= (1.0 + var_h) * hm[1]

      # KKT stationarity
      Lagr[mc=1:nv, i=1:nfe],
          d[mc] + w*v[mc,i]*vs[mc] +
          alpha_L[mc,i] + alpha_U[mc,i] +
          sum(SELECT_UPTAKE[mc,k] * alpha_upt[k,i] for k in 1:n_up) +
          sum(SELECT_PRODUCT[mc,k] * alpha_prod[k,i] for k in 1:n_prod) +
          sum(S[k,mc] * lambda_[k,i] for k in 1:nm) == 0

      # Flux smoothing between FE (paper-inspired)
      flux_smooth_pos[mc=1:nv, i=2:nfe],  v[mc,i] - v[mc,i-1] <= eps_flux
      flux_smooth_neg[mc=1:nv, i=2:nfe],  v[mc,i-1] - v[mc,i] <= eps_flux
  end)

  # --- Nonlinear constraints: ODEs + coupling + complementarity
  @NLconstraints(m, begin
      # ODEs (6 states)
      m1[i=1:nfe, j=1:ncp], cdot[1,i,j] ==
          ((v[obj,i]*vs[obj] - Kd_j[i,j]) * cX[i,j]) / cs[1]

      m2[i=1:nfe, j=1:ncp], cdot[2,i,j] ==
          - MW_N * (sum(IS_NIT[k] *
                        (-v[UPTAKE_IDXS[k],i] * vs[UPTAKE_IDXS[k]]) *
                        N_atoms_vec[UPTAKE_IDXS[k]]
                    for k in 1:n_up)) * cX[i,j] / cs[2] + injection_rate[i,j] / cs[2]

      m3[i=1:nfe, j=1:ncp], cdot[3,i,j] ==
          - MW_GLU * (-v[glu,i]*vs[glu]) * cX[i,j] / cs[3]

      m4[i=1:nfe, j=1:ncp], cdot[4,i,j] ==
          - MW_FRU * (-v[fru,i]*vs[fru]) * cX[i,j] / cs[4]

      m5[i=1:nfe, j=1:ncp], cdot[5,i,j] ==
          MW_ETH * v[eth,i]*vs[eth] * cX[i,j] / cs[5]

      # O2 por balance dinamico (sin cota cinetica explicita)
      m6[i=1:nfe, j=1:ncp], cdot[6,i,j] ==
          - MW_O2 * (-v[o2,i]*vs[o2]) * cX[i,j] / cs[6]

      # Production bounds (j=ncp)
      v_UB_product[k=1:n_prod, i=1:nfe],
          v[PRODUCT_IDXS[k],i]*vs[PRODUCT_IDXS[k]] - L_product_cap[k,i] <= 0

      # Uptake bounds (j=ncp)
      v_LB_uptake[k=1:n_up, i=1:nfe],
          -v[UPTAKE_IDXS[k],i]*vs[UPTAKE_IDXS[k]] - L_uptake[k,i,ncp] <= 0

      # Complementarity
      FO1[mc=1:nv, i=1:nfe],
          FO_L[mc,i] == (v[mc,i]*vs[mc] - vlb[mc]) * alpha_L[mc,i]

      FO2[mc=1:nv, i=1:nfe],
          FO_U[mc,i] == (v[mc,i]*vs[mc] - vub[mc]) * alpha_U[mc,i]

      FO3[k=1:n_up, i=1:nfe],
          FO_upt[k,i] == (-v[UPTAKE_IDXS[k],i]*vs[UPTAKE_IDXS[k]] -
                           L_uptake[k,i,ncp]) * alpha_upt[k,i]

      FO4[k=1:n_prod, i=1:nfe],
          FO_prod[k,i] == (v[PRODUCT_IDXS[k],i]*vs[PRODUCT_IDXS[k]] -
                           L_product_cap[k,i]) * alpha_prod[k,i]
  end)

  JuMP.optimize!(m)
  term = termination_status(m)
  prim = primal_status(m)
  dual = dual_status(m)
  println("\nStatus = ", term)
  println("Primal = ", prim)

  cStar   = JuMP.value.(c[:,:,:])  .* cs
  vStar   = JuMP.value.(v[:,:])    .* vs
  cdStar  = JuMP.value.(cdot[:,:,:])
  lStar   = JuMP.value.(lambda_[:,:])
  alLStar = JuMP.value.(alpha_L[:,:])
  alUStar = JuMP.value.(alpha_U[:,:])
  agStar  = JuMP.value.(alpha_upt[:,:])
  apStar  = JuMP.value.(alpha_prod[:,:])
  hStar   = JuMP.value.(hv[:])
  L_upt_star = JuMP.value.(L_uptake[:,:,ncp])
  L_prod_star = JuMP.value.(L_product[:,:,ncp])
    L_prod_cap_star = JuMP.value.(L_product_cap[:,:])
  max_upt_viol = maximum(-vStar[UPTAKE_IDXS[k],i] - L_upt_star[k,i] for k in 1:n_up, i in 1:nfe)
    max_prod_viol = maximum(vStar[PRODUCT_IDXS[k],i] - L_prod_cap_star[k,i] for k in 1:n_prod, i in 1:nfe)
  diagStar = (
      vx = JuMP.value.(vx[:,ncp]),
      vg = JuMP.value.(vg[:,ncp]),
      vf = JuMP.value.(vf[:,ncp]),
      vn = JuMP.value.(vn[:,ncp]),
      ve = JuMP.value.(ve[:,ncp]),
      L_upt = L_upt_star,
      L_prod = L_prod_star,
    L_prod_cap = L_prod_cap_star,
    apply_product_caps = apply_product_caps,
      solver = (term=term, primal=prim, dual=dual),
      o2_abs_max = maximum(abs.(vStar[o2,:])),
      max_uptake_violation = max_upt_viol,
      max_product_violation = max_prod_viol,
  )

  return cStar, vStar, cdStar, lStar, alLStar, alUStar, agStar, apStar, hStar, diagStar
end
