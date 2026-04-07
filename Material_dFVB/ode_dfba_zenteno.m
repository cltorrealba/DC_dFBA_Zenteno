%% =========================================================================
%  ode_dfba_zenteno.m
%  =========================================================================
%  ODE system for dFBA simulation of wine fermentation (S. cerevisiae)
%  Based on Zenteno et al. kinetic model + genome-scale yeast-GEM
%
%  USAGE:
%    dy = ode_dfba_zenteno(t, y, params, v_fba)
%
%  INPUTS:
%    t      - current time [h]
%    y      - state vector (n_state x 1)
%    params - struct with kinetic parameters (see load_params.m)
%    v_fba  - struct with FBA fluxes at current state:
%               .v_obj, .v_glu, .v_fru, .v_eth, .v_o2, .v_atpm, .v_prot
%               .vn_eff     (effective N uptake, gN/gDW/h)
%               .v_aa(j)    (AA fluxes, mmol/gDW/h, negative = uptake)
%               .v_aroma(j) (aroma fluxes, mmol/gDW/h, positive = production)
%               .mode       (1 = BIOMASS, 2 = TURNOVER_ATPM)
%
%  OUTPUTS:
%    dy     - state derivative vector (n_state x 1)
%
%  STATE VECTOR y (indices):
%    1  X       Biomass                 [gDW/L]
%    2  N_free  Free nitrogen (NH4+)    [gN/L]
%    3  G       Glucose                 [g/L]
%    4  F       Fructose                [g/L]
%    5  E       Ethanol                 [g/L]
%    6  O2      Dissolved oxygen        [g/L]
%    7  Prot    Intracellular protein   [g/L]
%    8  Carb    Intracellular carbohydrate [g/L]
%    9..8+n_aa     AA_i  Amino acids in medium [mmol/L]
%    9+n_aa..end   Aroma_j  Higher alcohols    [g/L]
%  =========================================================================

function dy = ode_dfba_zenteno(t, y, params, v_fba)

    n_state = length(y);
    dy = zeros(n_state, 1);
    
    %% Unpack state
    cX    = max(y(1), 0);
    cN    = max(y(2), 0);
    cG    = max(y(3), 0);
    cF    = max(y(4), 0);
    cE    = max(y(5), 0);
    cO2   = max(y(6), 0);
    cProt = max(y(7), 0);
    cCarb = max(y(8), 0);
    
    %% Unpack FBA fluxes
    v_obj  = v_fba.v_obj;
    v_glu  = v_fba.v_glu;   % negative = uptake
    v_fru  = v_fba.v_fru;   % negative = uptake
    v_eth  = v_fba.v_eth;   % positive = production
    v_o2   = v_fba.v_o2;    % negative = uptake
    v_prot = v_fba.v_prot;
    vn_eff = v_fba.vn_eff;  % gN/gDW/h (positive = uptake)
    mode   = v_fba.mode;    % 1=BIOMASS, 2=TURNOVER
    
    %% Extract parameters
    MW_GLU  = params.MW_GLU;   % 0.180156 g/mmol
    MW_FRU  = params.MW_FRU;   % 0.180156 g/mmol
    MW_ETH  = params.MW_ETH;   % 0.046070 g/mmol
    MW_O2   = params.MW_O2;    % 0.032    g/mmol
    K_DEATH = params.K_DEATH;  % 0.005    1/h
    TURNOVER_LAMBDA = params.TURNOVER_LAMBDA; % 0.03 1/h
    PROT_CONTENT_0  = params.PROT_CONTENT_0;  % 0.46 g/gDW
    CARB_CONTENT_0  = params.CARB_CONTENT_0;  % 0.37 g/gDW
    XA_FRACTION     = params.XA_FRACTION;     % 1.0
    
    n_aa    = params.n_aa;
    n_aroma = params.n_aroma;
    
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
    dy(8) = max(0, v_obj) * cX * CARB_CONTENT_0 ...
            - cCarb * K_DEATH;
    
    % (9..8+n_aa) Amino acids [mmol/L]
    for j = 1:n_aa
        idx = 8 + j;
        if mode == 1  % BIOMASS
            dy(idx) = cX * v_fba.v_aa(j);  % v_aa negative = uptake
        else
            dy(idx) = 0;
        end
    end
    
    % (8+n_aa+1..end) Aromas [g/L]
    for j = 1:n_aroma
        idx = 8 + n_aa + j;
        mw_j = params.aroma_MW(j);  % g/mmol
        dy(idx) = mw_j * max(0, v_fba.v_aroma(j)) * cX;
    end
    
end
