%% =========================================================================
%  kinetic_limits_zenteno.m
%  =========================================================================
%  Compute kinetic upper bounds for FBA fluxes based on Zenteno model.
%  These bounds are used to constrain the FBA problem at each time step.
%
%  USAGE:
%    limits = kinetic_limits_zenteno(cX, cN, cG, cF, cE, cO2, params)
%
%  INPUTS:
%    cX   - biomass concentration [gDW/L]
%    cN   - free nitrogen [gN/L]
%    cG   - glucose [g/L]
%    cF   - fructose [g/L]
%    cE   - ethanol [g/L]
%    cO2  - dissolved oxygen [g/L]
%    params - struct with kinetic parameters
%
%  OUTPUTS:
%    limits - struct with fields:
%      .lim_glu  - max glucose uptake [mmol/gDW/h]
%      .lim_fru  - max fructose uptake [mmol/gDW/h]
%      .lim_eth  - max ethanol production [mmol/gDW/h]
%      .lim_obj  - max growth rate [1/h]
%      .lim_n    - max N uptake per source (vector) [mmol/gDW/h]
%  =========================================================================

function limits = kinetic_limits_zenteno(cX, cN, cG, cF, cE, cO2, params)
    
    EPS = params.EPS;      % 1e-9
    T   = params.T_val;    % 293.15 K
    R   = params.R_GAS;    % 8.314 J/mol/K
    
    % Temperature corrections (Arrhenius-type)
    mu_T  = exp(59453.0 * (T - 300.0)  / (300.0  * R * T));
    Kg_T  = exp(46055.0 * (T - 293.15) / (293.15 * R * T));
    b_T   = exp(11000.0 * (T - 296.15) / (296.15 * R * T));
    
    MRATE_0 = params.MRATE_0;
    mrate = MRATE_0 * exp(37681.0 * (T - 293.30) / (293.30 * R * T));
    
    % Zenteno kinetic equations
    % Growth rate (Monod on N)
    mu = params.MU0_nom * mu_T * (cN / (cN + params.Kn0_nom * Kg_T + EPS));
    
    % Glucose uptake (Monod + ethanol inhibition)
    betaG = params.betaG0_nom * b_T ...
            * (cG / (cG + params.Kg0_nom * Kg_T + EPS)) ...
            * (params.Kie0_nom * Kg_T / (cE + params.Kie0_nom * Kg_T + EPS));
    
    % Fructose uptake (Monod + glucose inhibition + ethanol inhibition)
    betaF = params.betaF0_nom * b_T ...
            * (cF / (cF + params.Kf0_nom * Kg_T + EPS)) ...
            * (params.Kig0_nom * Kg_T / (cG + params.Kig0_nom * Kg_T + EPS)) ...
            * (params.Kie0_nom * Kg_T / (cE + params.Kie0_nom * Kg_T + EPS));
    
    % Derived flux limits
    vx = mu;
    vg = mu / params.YXG_nom + betaG / params.YEG_nom ...
         + mrate * (cG / (cG + cF + EPS));
    vf = mu / params.YXF_nom + betaF / params.YEF_nom ...
         + mrate * (cF / (cG + cF + EPS));
    vn = mu / params.YXN_nom;
    ve = (betaG + betaF) / params.MW_ETH;
    
    % Convert to mmol/gDW/h (flux units)
    limits.lim_glu = max(0, vg / params.MW_GLU);
    limits.lim_fru = max(0, vf / params.MW_FRU);
    limits.lim_eth = max(0, ve);
    limits.lim_obj = max(0, vx);
    
    % N uptake per source (distributed by N_frac_share)
    vn_total = max(0, vn);
    limits.lim_n = vn_total * params.N_frac_share;  % vector
    
end
