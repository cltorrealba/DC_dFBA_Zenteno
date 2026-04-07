%% =========================================================================
%  compute_full_gam.m
%  =========================================================================
%  Calculate variable Growth-Associated Maintenance (GAM) based on
%  dynamic protein, RNA, and carbohydrate composition.
%  Based on Sanchez et al. (2017) changeGAM approach.
%
%  USAGE:
%    fullGAM = compute_full_gam(P, R, C, Pbase, Rbase, Cbase, params)
%
%  INPUTS:
%    P      - current protein fraction [g protein / gDW]
%    R      - current RNA fraction [g RNA / gDW] (usually fixed = 0.06)
%    C      - current carbohydrate fraction [g carb / gDW]
%    Pbase  - reference protein fraction [g/gDW] (from yeast-GEM biomass)
%    Rbase  - reference RNA fraction [g/gDW]
%    Cbase  - reference carbohydrate fraction [g/gDW]
%    params - struct with GAM coefficients
%
%  OUTPUTS:
%    fullGAM - Growth-associated maintenance [mmol ATP / gDW]
%
%  FORMULA:
%    fullGAM = GAM_BASE + GAM_COEFF_P * (P/Pbase) 
%                       + GAM_COEFF_R * (R/Rbase) 
%                       + GAM_COEFF_C * max(0, (Cbase+Pbase-P-R)/Cbase)
%  =========================================================================

function fullGAM = compute_full_gam(P, R, C, Pbase, Rbase, Cbase, params)

    EPS = 1e-9;
    
    Pfactor = P / max(Pbase, EPS);
    Rfactor = R / max(Rbase, EPS);
    Cfactor = (Cbase + Pbase - P - R) / max(Cbase, EPS);
    Cfactor = max(0, Cfactor);
    
    fullGAM = params.GAM_BASE ...
              + params.GAM_COEFF_P * Pfactor ...
              + params.GAM_COEFF_R * Rfactor ...
              + params.GAM_COEFF_C * Cfactor;

end
