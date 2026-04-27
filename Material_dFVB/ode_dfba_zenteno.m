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
%               .v_aroma(j) (fallback aroma fluxes, mmol/gDW/h, positive)
%               .mode       (1 = BIOMASS, 2 = TURNOVER_ATPM)
%             Aromas with Henriques product_k are integrated directly as
%               vP_i = k_i * fermentation_hexose_Zenteno(T) * f_AA_i
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
    MW_N    = get_param_value(params, {'MW_N'}, 0.014007); % g/mmol N
    EPS     = get_param_value(params, {'EPS'}, 1e-9);
    K_DEATH = params.K_DEATH;  % 0.005    1/h
    TURNOVER_LAMBDA = params.TURNOVER_LAMBDA; % 0.03 1/h
    PROT_CONTENT_0  = params.PROT_CONTENT_0;  % 0.46 g/gDW
    CARB_CONTENT_0  = params.CARB_CONTENT_0;  % 0.37 g/gDW
    XA_FRACTION     = params.XA_FRACTION;     % 1.0
    
    n_aa    = params.n_aa;
    n_aroma = params.n_aroma;
    aroma_names = get_aroma_names(params, n_aroma);
    aroma_MW = get_aroma_mw(params, n_aroma);
    aroma_product_k = get_aroma_product_k(params, aroma_names);
    use_henriques_aroma = get_param_bool(params, {'enable_henriques_aroma_ode'}, true);
    fermentation_hexose = zenteno_fermentation_hexose( ...
        t, cG, cF, cE, params, MW_GLU, MW_FRU, EPS);
    aroma_product_factors = get_aroma_product_factors( ...
        y, n_aa, aroma_names, params, MW_N, EPS);
    
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
        mw_j = aroma_MW(j);  % g/mmol
        use_direct_rate = use_henriques_aroma && aroma_product_k(j) > 0;
        if use_direct_rate
            % Notebook update: vP_i = k_i * fermentation_hexose * f_AA_i.
            v_aroma_j = aroma_product_k(j) ...
                        * max(0, fermentation_hexose) ...
                        * aroma_product_factors(j);
            dy(idx) = mw_j * max(0, v_aroma_j) * Xa;
        else
            dy(idx) = mw_j * max(0, get_vector_field(v_fba, 'v_aroma', j, 0)) * cX;
        end
    end
    
end

function fermentation_hexose = zenteno_fermentation_hexose(t, cG, cF, cE, params, MW_GLU, MW_FRU, EPS)
    T = temperature_at_time(t, params);
    R = get_param_value(params, {'R_GAS'}, 8.314);

    Kg_T = exp(46055.0 * (T - 293.15) / (293.15 * R * T));
    b_T  = exp(11000.0 * (T - 296.15) / (296.15 * R * T));

    betaG0 = get_param_value(params, {'betaG0_nom', 'betaG'}, 0.225);
    betaF0 = get_param_value(params, {'betaF0_nom', 'betaF'}, 0.225);
    Kg0    = get_param_value(params, {'Kg0_nom', 'Kg'}, 7.5);
    Kf0    = get_param_value(params, {'Kf0_nom', 'Kf'}, 7.5);
    Kig0   = get_param_value(params, {'Kig0_nom', 'Kig'}, 55.0);
    Kie0   = get_param_value(params, {'Kie0_nom', 'Kie'}, 40.0);
    YEG    = get_param_value(params, {'YEG_nom', 'YEG'}, 0.49);
    YEF    = get_param_value(params, {'YEF_nom', 'YEF'}, 0.49);

    betaG = betaG0 * b_T ...
            * (cG / (cG + Kg0 * Kg_T + EPS)) ...
            * (Kie0 * Kg_T / (cE + Kie0 * Kg_T + EPS));

    betaF = betaF0 * b_T ...
            * (cF / (cF + Kf0 * Kg_T + EPS)) ...
            * (Kig0 * Kg_T / (cG + Kig0 * Kg_T + EPS)) ...
            * (Kie0 * Kg_T / (cE + Kie0 * Kg_T + EPS));

    fermentation_hexose = betaG / max(YEG * MW_GLU, EPS) ...
                        + betaF / max(YEF * MW_FRU, EPS);
end

function T = temperature_at_time(t, params)
    profile = get_param_value(params, {'temperature_profile_df', 'temperature_profile'}, []);
    if istable(profile) ...
            && all(ismember({'t_h', 'T_K'}, profile.Properties.VariableNames)) ...
            && height(profile) > 0
        t_arr = profile.t_h;
        T_arr = profile.T_K;
        [t_arr, order] = sort(t_arr);
        T_arr = T_arr(order);
        mode = get_param_value(params, {'temperature_profile_mode'}, 'step');
        if strcmpi(char(mode), 'linear')
            T = interp1(t_arr, T_arr, t, 'linear', 'extrap');
        else
            idx = find(t_arr <= t, 1, 'last');
            if isempty(idx)
                idx = 1;
            end
            T = T_arr(idx);
        end
        return;
    end

    has_dynamic_steps = has_param_value(params, 'T_base_K') ...
                        && has_param_value(params, 'T_steps_h') ...
                        && has_param_value(params, 'T_deltas_K');
    if has_dynamic_steps
        T = get_param_value(params, {'T_base_K'}, 298.15);
        steps = get_param_value(params, {'T_steps_h'}, []);
        deltas = get_param_value(params, {'T_deltas_K'}, []);
        steep = get_param_value(params, {'T_steep'}, 0.5);
        for k = 1:min(numel(steps), numel(deltas))
            sigmoid = 1.0 / (1.0 + exp(-steep * (t - steps(k))));
            T = T + deltas(k) * sigmoid;
        end
    else
        T = get_param_value(params, {'T_val', 'T_base_K', 'T'}, 293.15);
    end
end

function factors = get_aroma_product_factors(y, n_aa, aroma_names, params, MW_N, EPS)
    n_aroma = numel(aroma_names);
    factors = ones(1, n_aroma);
    if ~get_param_bool(params, {'enable_aroma_aa_modulation'}, true)
        return;
    end

    aa_names = get_aa_names(params, n_aa);
    aa_factors = struct();
    km = max(get_param_value(params, {'aroma_aa_K_gN_L'}, 0.002), EPS);
    floor_val = min(max(get_param_value(params, {'aroma_aa_floor'}, 0.05), 0.0), 1.0);

    for j = 1:n_aa
        aa_key = matlab.lang.makeValidName(aa_names{j});
        c_aa_mmol_L = max(0, y(8 + j));
        c_aa_gN_L = c_aa_mmol_L * get_aa_n_atoms(params, aa_names{j}) * MW_N;
        sat = c_aa_gN_L / (c_aa_gN_L + km + EPS);
        aa_factors.(aa_key) = floor_val + (1.0 - floor_val) * sat;
    end

    for j = 1:n_aroma
        deps = get_product_dependencies(params, aroma_names{j});
        factors(j) = combine_dependency_factors(deps, aa_factors);
    end
end

function value = combine_dependency_factors(deps, aa_factors)
    value = 1.0;
    weight_sum = 0.0;
    dep_names = fieldnames(deps);
    for k = 1:numel(dep_names)
        aa_key = matlab.lang.makeValidName(dep_names{k});
        weight = max(0.0, deps.(dep_names{k}));
        if weight <= 0.0
            continue;
        end
        if isfield(aa_factors, aa_key)
            value = value * aa_factors.(aa_key) ^ weight;
            weight_sum = weight_sum + weight;
        end
    end
    if weight_sum <= 0.0
        value = 1.0;
    end
end

function deps = get_product_dependencies(params, product_name)
    deps = struct();
    dep_map = get_param_value(params, {'aroma_product_aa_dependency'}, []);
    product_field = matlab.lang.makeValidName(product_name);
    if isstruct(dep_map) && isfield(dep_map, product_field)
        deps = dep_map.(product_field);
        return;
    end

    switch canonical_aroma_name(product_name)
        case {'isoamyl_acetate', 'isoamyl_alcohol'}
            deps.leu = 1.0;
        case {'isobutyl_acetate'}
            deps.val = 1.0;
        case {'methionol'}
            deps.met = 1.0;
        case {'phenylethanol'}
            deps.phe = 1.0;
    end
end

function product_k = get_aroma_product_k(params, aroma_names)
    n_aroma = numel(aroma_names);
    product_k = zeros(1, n_aroma);
    configured = get_param_value(params, {'aroma_product_k', 'product_k'}, []);

    if isnumeric(configured) && numel(configured) >= n_aroma
        product_k = configured(1:n_aroma);
        product_k = reshape(product_k, 1, []);
        return;
    end

    for j = 1:n_aroma
        product_name = canonical_aroma_name(aroma_names{j});
        product_k(j) = product_k_for_name(configured, product_name);
    end
end

function aroma_MW = get_aroma_mw(params, n_aroma)
    aroma_MW = get_param_value(params, {'aroma_MW'}, []);
    if isnumeric(aroma_MW) && numel(aroma_MW) >= n_aroma
        aroma_MW = reshape(aroma_MW(1:n_aroma), 1, []);
        return;
    end

    aroma_MW = zeros(1, n_aroma);
    if isfield(params, 'aroma') && istable(params.aroma) ...
            && ismember('MW_g_per_mmol', params.aroma.Properties.VariableNames)
        for j = 1:min(n_aroma, height(params.aroma))
            aroma_MW(j) = params.aroma.MW_g_per_mmol(j);
        end
    end
end

function value = product_k_for_name(configured, product_name)
    value = NaN;
    field = matlab.lang.makeValidName(product_name);
    if isstruct(configured) && isfield(configured, field)
        value = configured.(field);
    end

    if isnan(value)
        switch product_name
            case 'isoamyl_acetate'
                value = 2.0e-5;
            case 'ethyl_acetate'
                value = 2.0e-4;
            otherwise
                value = 0.0;
        end
    end
end

function names = get_aroma_names(params, n_aroma)
    defaults = {'phenylethanol', 'isoamyl_acetate', 'isobutanol', ...
                'methionol', 'tyrosol', 'ethyl_acetate'};
    names = default_names(defaults, n_aroma, 'aroma');

    if isfield(params, 'aroma') && istable(params.aroma)
        T = params.aroma;
        for j = 1:min(n_aroma, height(T))
            raw_key = '';
            raw_rid = '';
            raw_rxn = '';
            if ismember('aroma_key', T.Properties.VariableNames)
                raw_key = table_value_as_char(T, 'aroma_key', j);
            end
            if ismember('reaction_id', T.Properties.VariableNames)
                raw_rid = table_value_as_char(T, 'reaction_id', j);
            end
            if ismember('reaction_name', T.Properties.VariableNames)
                raw_rxn = table_value_as_char(T, 'reaction_name', j);
            end
            raw = sprintf('%s %s %s', raw_key, raw_rid, raw_rxn);
            names{j} = canonical_aroma_name(raw);
        end
    elseif isfield(params, 'aroma_keys')
        names = normalize_names(params.aroma_keys, names, n_aroma);
    end
end

function names = get_aa_names(params, n_aa)
    defaults = {'phe', 'leu', 'val', 'met', 'tyr'};
    names = default_names(defaults, n_aa, 'aa');

    if isfield(params, 'aa') && istable(params.aa) ...
            && ismember('aa_key', params.aa.Properties.VariableNames)
        for j = 1:min(n_aa, height(params.aa))
            names{j} = canonical_aa_name(table_value_as_char(params.aa, 'aa_key', j));
        end
    elseif isfield(params, 'aa_keys')
        names = normalize_names(params.aa_keys, names, n_aa);
        for j = 1:n_aa
            names{j} = canonical_aa_name(names{j});
        end
    end
end

function n_atoms = get_aa_n_atoms(params, aa_name)
    field = matlab.lang.makeValidName(aa_name);
    n_atoms = 1.0;
    configured = get_param_value(params, {'aa_N_atoms', 'N_atoms'}, []);
    if isstruct(configured) && isfield(configured, field)
        n_atoms = configured.(field);
    end
end

function names = default_names(defaults, n, prefix)
    names = cell(1, n);
    for j = 1:n
        if j <= numel(defaults)
            names{j} = defaults{j};
        else
            names{j} = sprintf('%s_%d', prefix, j);
        end
    end
end

function names = normalize_names(raw_names, fallback, n)
    names = fallback;
    if isstring(raw_names) || iscell(raw_names)
        for j = 1:min(n, numel(raw_names))
            names{j} = value_as_char(raw_names(j));
        end
    end
end

function name = canonical_aroma_name(raw_name)
    s = lower(regexprep(value_as_char(raw_name), '[^a-z0-9]+', '_'));
    if contains(s, 'r_1765') || (contains(s, 'ethyl') && contains(s, 'acetate'))
        name = 'ethyl_acetate';
    elseif contains(s, 'r_1862') || contains(s, 'isoamyl')
        name = 'isoamyl_acetate';
    elseif contains(s, 'pea') || contains(s, 'phenylethanol') || contains(s, 'phenyl_ethanol')
        name = 'phenylethanol';
    elseif contains(s, 'isobutyl') && contains(s, 'acetate')
        name = 'isobutyl_acetate';
    elseif contains(s, 'isobutanol')
        name = 'isobutanol';
    elseif contains(s, 'methionol')
        name = 'methionol';
    elseif contains(s, 'tyrosol')
        name = 'tyrosol';
    else
        name = s;
    end
end

function name = canonical_aa_name(raw_name)
    s = lower(regexprep(value_as_char(raw_name), '[^a-z0-9]+', '_'));
    if contains(s, 'phenylalanine') || strcmp(s, 'phe')
        name = 'phe';
    elseif contains(s, 'leucine') || strcmp(s, 'leu') || contains(s, 'isoleucine')
        name = 'leu';
    elseif contains(s, 'valine') || strcmp(s, 'val')
        name = 'val';
    elseif contains(s, 'methionine') || strcmp(s, 'met')
        name = 'met';
    elseif contains(s, 'tyrosine') || strcmp(s, 'tyr')
        name = 'tyr';
    else
        name = s;
    end
end

function value = get_vector_field(s, field_name, idx, default_value)
    value = default_value;
    if isstruct(s) && isfield(s, field_name)
        vec = s.(field_name);
        if numel(vec) >= idx
            value = vec(idx);
        end
    end
end

function value = get_param_value(params, names, default_value)
    value = default_value;
    if ischar(names) || isstring(names)
        names = cellstr(names);
    end

    for k = 1:numel(names)
        name = char(names{k});
        if isstruct(params) && isfield(params, name)
            value = params.(name);
            return;
        end
    end

    nested_blocks = {'kinetic', 'turnover'};
    for b = 1:numel(nested_blocks)
        block = nested_blocks{b};
        if isstruct(params) && isfield(params, block) && isstruct(params.(block))
            for k = 1:numel(names)
                name = char(names{k});
                if isfield(params.(block), name)
                    value = params.(block).(name);
                    return;
                end
            end
        end
    end
end

function tf = has_param_value(params, name)
    sentinel = '__missing__';
    tf = ~isequal(get_param_value(params, {name}, sentinel), sentinel);
end

function value = get_param_bool(params, names, default_value)
    value = get_param_value(params, names, default_value);
    if isnumeric(value) || islogical(value)
        value = logical(value);
    elseif ischar(value) || isstring(value)
        value = any(strcmpi(char(value), {'true', '1', 'yes', 'on'}));
    else
        value = default_value;
    end
end

function value = table_value_as_char(T, column_name, idx)
    column = T.(column_name);
    if iscell(column)
        value = value_as_char(column{idx});
    elseif isstring(column)
        value = char(column(idx));
    elseif iscategorical(column)
        value = char(column(idx));
    else
        value = value_as_char(column(idx));
    end
end

function value = value_as_char(raw)
    if iscell(raw)
        value = value_as_char(raw{1});
    elseif isstring(raw)
        value = char(raw(1));
    elseif ischar(raw)
        value = raw;
    elseif isnumeric(raw)
        value = num2str(raw(1));
    else
        value = '';
    end
end
