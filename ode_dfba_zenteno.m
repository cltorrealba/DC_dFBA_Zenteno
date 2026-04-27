function dy = ode_dfba_zenteno(t, y, model, params, adj_params)
    
    persistent sol_factible
    if t == 0
        sol_factible = struct();
    elseif strcmp(sol_factible.status,'OPTIMAL')
        model.vbasis = sol_factible.vbasis;
        model.cbasis = sol_factible.cbasis;
    end
    
    gurobiparams.OutputFlag  = 0;
    gurobiparams.LPWarmStart = 1;

    % Unpack state
    cX    = y(1);
    cN    = y(2);
    cG    = y(3);
    cF    = y(4);
    cE    = y(5);
    cO2   = y(6);
    cProt = y(7);
    cCarb = y(8);
    AA    = y(9:end);

    % Extract parameters
    MW_GLU             = 0.180156; % g/mmol
    MW_FRU             = 0.180156; % g/mmol
    MW_ETH             = 0.046070; % g/mmol
    MW_O2              = 0.032;    % g/mmol
    MW_N               = 0.014007; % g/mmol
    K_DEATH            = 0.005;    % 1/h
    TURNOVER_LAMBDA    = 0.03;     % 1/h
    PROT_CONTENT_0     = 0.46;     % g/gDW
    CARB_CONTENT_0     = 0.37;     % g/gDW
    XA_FRACTION        = 1.0;   
    n_aa               = 5;
    n_aroma            = 6;
    EPS                = 1e-9;      
    T                  = 293.15;   % K
    R                  = 8.314;    % J/mol/K
    TURNOVER_THRESHOLD = 0.001; % g/L
    ATP_LB_TURNOVER    = 0.7;
    ATP_UB_TURNOVER    = 1e3;

    % Adjustable parameters
    MU0_nom    = adj_params(1);
    YXN_nom    = adj_params(2);
    YXG_nom    = adj_params(3);
    YXF_nom    = adj_params(4);
    YEG_nom    = adj_params(5);
    YEF_nom    = adj_params(6);
    Kn0_nom    = adj_params(7);
    Kg0_nom    = adj_params(8);
    Kf0_nom    = adj_params(9);
    Kig0_nom   = adj_params(10);
    Kie0_nom   = adj_params(11);
    Kd0_nom    = adj_params(12);
    betaG0_nom = adj_params(13);
    betaF0_nom = adj_params(14);

    % Temperature corrections (Arrhenius-type)
    mu_T  = exp(59453.0 * (T - 300.0)  / (300.0  * R * T));
    Kg_T  = exp(46055.0 * (T - 293.15) / (293.15 * R * T));
    b_T   = exp(11000.0 * (T - 296.15) / (296.15 * R * T));
    
    MRATE_0 = 0.01;
    mrate = MRATE_0 * exp(37681.0 * (T - 293.30) / (293.30 * R * T));
    
    % Zenteno kinetic equations
    % Growth rate (Monod on N)
    mu = MU0_nom * mu_T * (cN / (cN + Kn0_nom * Kg_T + EPS));
    
    % Glucose uptake (Monod + ethanol inhibition)
    betaG = betaG0_nom * b_T ...
            * (cG / (cG + Kg0_nom * Kg_T + EPS)) ...
            * (Kie0_nom * Kg_T / (cE + Kie0_nom * Kg_T + EPS));
    
    % Fructose uptake (Monod + glucose inhibition + ethanol inhibition)
    betaF = betaF0_nom * b_T ...
            * (cF / (cF + Kf0_nom * Kg_T + EPS)) ...
            * (Kig0_nom * Kg_T / (cG + Kig0_nom * Kg_T + EPS)) ...
            * (Kie0_nom * Kg_T / (cE + Kie0_nom * Kg_T + EPS));
    
    % Derived flux limits
    vx = mu;
    vg = mu / YXG_nom + betaG / YEG_nom + mrate * (cG / (cG + cF + EPS));
    vf = mu / YXF_nom + betaF / YEF_nom + mrate * (cF / (cG + cF + EPS));
    vn = mu / YXN_nom;
    ve = (betaG + betaF) / MW_ETH;
       
    % N uptake per source (distributed by N_frac_share)
    vn_total = max(0, vn);

    % FBA

    % Turnover
    N_total = cN + sum(AA) * MW_N;
    if N_total > TURNOVER_THRESHOLD
        mode = 1;
        model.c = params.objective_turnover_1;
        model.lb(params.atp_rxn) = ATP_LB_TURNOVER;
        model.ub(params.atp_rxn) = ATP_UB_TURNOVER;
    else
        mode = 0;
        model.lb(params.growth_rxn) = 0;
        model.ub(params.growth_rxn) = 0;
        model.c = params.objective_turnover_0;
    end

    model.lb(params.glc_rxn)    = max(0, vg / MW_GLU);
    model.lb(params.fru_rxn)    = max(0, vf / MW_FRU);
    model.lb(params.eth_rxn)    = max(0, ve);
    model.lb(params.growth_rxn) = max(0, vx);
    model.lb(params.idx_n)      = vn_total * params.N_frac_share;
    
    sol = gurobi(model, gurobiparams);
    sol_factible = sol;

    n_state = length(y);
    dy = zeros(n_state, 1);

    if strcmp(sol.status, 'OPTIMAL')        
        
        %% Unpack FBA fluxes
        v_obj    = sol.x(params.growth_rxn);
        v_glu    = sol.x(params.glc_rxn);   % negative = uptake
        v_fru    = sol.x(params.fru_rxn);   % negative = uptake
        v_eth    = sol.x(params.eth_rxn);   % positive = production
        v_o2     = sol.x(params.o2_rxn);    % negative = uptake
        v_prot   = sol.x(params.prot_rxn);
        vn_eff   = sol.x(params.vn_eff_rxn);  % gN/gDW/h (positive = uptake)
        v_aa     = sol.x(params.aa_rxn);
        v_aromas = sol.x(params.aromas_rxn);

        %% === ODE equations ===
        
        % (1) Biomass [gDW/L]
        if mode == 1  % BIOMASS
            dy(1) = v_obj * cX;
        else          % TURNOVER
            dy(1) = 0;
        end
        
        % (2) Free nitrogen [gN/L]
        if mode == 1
            dy(2) = -vn_eff * cX;
        else
            dy(2) = 0;
        end
        
        % (3) Glucose [g/L]
        dy(3) = -MW_GLU * max(0, -v_glu) * cX;
        
        % (4) Fructose [g/L]
        dy(4) = -MW_FRU * max(0, -v_fru) * cX;
        
        % (5) Ethanol [g/L]
        dy(5) = MW_ETH * max(0, v_eth) * cX;
        
        % (6) O2 [g/L]
        dy(6) = -MW_O2 * max(0, -v_o2) * cX;
        
        % (7) Intracellular protein [g/L]
        Xa = XA_FRACTION * cX;
        dy(7) = max(0, v_obj) * cX * PROT_CONTENT_0 ...
                - cProt * K_DEATH ...
                + Xa * max(0, v_prot) ...
                - TURNOVER_LAMBDA * cProt;
        
        % (8) Intracellular carbohydrate [g/L]
        dy(8) = max(0, v_obj) * cX * CARB_CONTENT_0 - cCarb * K_DEATH;
        
        % (9..8+n_aa) Amino acids [mmol/L]
        for j = 1:n_aa
            idx = 8 + j;
            if mode == 1  % BIOMASS
                dy(idx) = cX * v_aa(j);  % v_aa negative = uptake
            else
                dy(idx) = 0;
            end
        end
        
        % (8+n_aa+1..end) Aromas [g/L]
        for j = 1:n_aroma
            idx = 8 + n_aa + j;
            mw_j = params.aroma_MW(j);  % g/mmol
            dy(idx) = mw_j * max(0, v_aromas(j)) * cX;
        end
    end
    
end
