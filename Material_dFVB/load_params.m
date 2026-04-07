%% =========================================================================
%  load_params.m
%  =========================================================================
%  Load all dFVB parameter files exported from Python (Notebook 1).
%  Returns a single struct 'params' with all kinetic, GAM, turnover,
%  objective, and model configuration data.
%
%  USAGE:
%    params = load_params(data_dir)
%
%  INPUT:
%    data_dir - path to the Material_dFVB/csv_exports/ folder
%
%  OUTPUT:
%    params - struct with fields:
%       .kinetic     - Zenteno kinetic parameters (MU0, YXG, etc.)
%       .gam         - GAM coefficients and base values
%       .turnover    - protein turnover / nitrogen recycling params
%       .nitrogen    - N_atoms_map, n_aa names, MW
%       .initial     - initial conditions (X0, G0, N0, F0, E0, O2_0, ...)
%       .objective   - objective function configuration per phase
%       .aroma       - aroma exchange reaction information
%       .aa          - amino acid exchange reaction names & indices
%       .ode_states  - ODE state definitions
%  =========================================================================

function params = load_params(data_dir)

    if nargin < 1
        data_dir = fullfile(pwd, 'csv_exports');
    end

    fprintf('Loading parameters from: %s\n', data_dir);

    % ---- Kinetic parameters ----
    T = readtable(fullfile(data_dir, 'kinetic_parameters.csv'));
    params.kinetic = struct();
    for i = 1:height(T)
        pname = strrep(T.parameter{i}, '.', '_');
        params.kinetic.(pname) = T.value(i);
    end

    % ---- GAM parameters ----
    T = readtable(fullfile(data_dir, 'gam_parameters.csv'));
    params.gam = struct();
    for i = 1:height(T)
        pname = T.parameter{i};
        params.gam.(pname) = T.value(i);
    end

    % ---- Turnover parameters ----
    T = readtable(fullfile(data_dir, 'turnover_parameters.csv'));
    params.turnover = struct();
    for i = 1:height(T)
        pname = T.parameter{i};
        params.turnover.(pname) = T.value(i);
    end

    % ---- Nitrogen configuration ----
    fpath = fullfile(data_dir, 'nitrogen_config.csv');
    if isfile(fpath)
        T = readtable(fpath);
        params.nitrogen = T;
    else
        params.nitrogen = [];
        warning('nitrogen_config.csv not found.');
    end

    % ---- Initial conditions ----
    T = readtable(fullfile(data_dir, 'initial_conditions.csv'));
    params.initial = struct();
    for i = 1:height(T)
        pname = T.state{i};
        params.initial.(pname) = T.value(i);
    end

    % ---- Objective configuration ----
    T = readtable(fullfile(data_dir, 'objective_config.csv'));
    params.objective = T;

    % ---- Aroma configuration ----
    fpath = fullfile(data_dir, 'aroma_config.csv');
    if isfile(fpath)
        T = readtable(fpath);
        params.aroma = T;
    else
        params.aroma = [];
        warning('aroma_config.csv not found.');
    end

    % ---- Amino acid configuration ----
    fpath = fullfile(data_dir, 'aa_precursor_config.csv');
    if isfile(fpath)
        T = readtable(fpath);
        params.aa = T;
    else
        params.aa = [];
        warning('aa_precursor_config.csv not found.');
    end

    % ---- ODE state definitions ----
    fpath = fullfile(data_dir, 'ode_state_definitions.csv');
    if isfile(fpath)
        T = readtable(fpath);
        params.ode_states = T;
    else
        params.ode_states = [];
        warning('ode_state_definitions.csv not found.');
    end

    fprintf('Parameters loaded successfully.\n');
end
