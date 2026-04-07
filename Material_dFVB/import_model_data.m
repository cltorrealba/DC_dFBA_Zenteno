%% =========================================================================
%  import_model_data.m
%  =========================================================================
%  Master import script: loads the stoichiometric matrix, bounds, reaction
%  & metabolite info, and all parameters into the MATLAB workspace.
%  
%  Run this first to set up everything for dFVB analysis.
%
%  OUTPUTS (in workspace):
%    S            - stoichiometric matrix [m x n]  (sparse)
%    rxn_info     - table with reaction id, name, subsystem, reversible
%    met_info     - table with metabolite id, name, compartment
%    lb           - lower bounds vector [n x 1]
%    ub           - upper bounds vector [n x 1]
%    rxn_indices  - table with reaction category indices
%    dyn_bounds   - table with dynamic bound rules
%    params       - struct with all parameters (via load_params)
%    bounds_log   - table with logged set_bounds calls (if available)
%  =========================================================================

%% Configuration
data_dir = fullfile(pwd, 'csv_exports');

fprintf('=== dFVB Model Import ===\n');
fprintf('Data directory: %s\n\n', data_dir);

%% 1. Stoichiometric matrix (sparse triplet format)
fprintf('[1/7] Loading stoichiometric matrix...\n');
T = readtable(fullfile(data_dir, 'S_matrix_sparse.csv'));
S = sparse(T.row + 1, T.col + 1, T.value);  % Python 0-indexed → MATLAB 1-indexed
fprintf('  S matrix: %d metabolites x %d reactions, %d nonzeros\n', ...
    size(S,1), size(S,2), nnz(S));

%% 2. Reaction info
fprintf('[2/7] Loading reaction information...\n');
rxn_info = readtable(fullfile(data_dir, 'reaction_info.csv'));
fprintf('  %d reactions loaded\n', height(rxn_info));

%% 3. Metabolite info
fprintf('[3/7] Loading metabolite information...\n');
met_info = readtable(fullfile(data_dir, 'metabolite_info.csv'));
fprintf('  %d metabolites loaded\n', height(met_info));

%% 4. Bounds
fprintf('[4/7] Loading bounds...\n');
lb = readmatrix(fullfile(data_dir, 'lb_vector.csv'));
ub = readmatrix(fullfile(data_dir, 'ub_vector.csv'));
fprintf('  lb: [%d x 1], ub: [%d x 1]\n', length(lb), length(ub));

% Sanity check
assert(length(lb) == size(S,2), 'lb length mismatch with S columns!');
assert(length(ub) == size(S,2), 'ub length mismatch with S columns!');

%% 5. Dynamic bounds
fprintf('[5/7] Loading dynamic bound rules...\n');
fpath = fullfile(data_dir, 'bounds_dynamic_rules.csv');
if isfile(fpath)
    dyn_bounds = readtable(fpath);
    fprintf('  %d dynamic bound rules loaded\n', height(dyn_bounds));
else
    dyn_bounds = [];
    fprintf('  (no dynamic_bounds.csv found)\n');
end

%% 6. Reaction indices
fprintf('[6/7] Loading reaction indices...\n');
rxn_indices = readtable(fullfile(data_dir, 'reaction_indices_of_interest.csv'));
fprintf('  Categories: ');
cats = unique(rxn_indices.category);
for k = 1:length(cats)
    n = sum(strcmp(rxn_indices.category, cats{k}));
    fprintf('%s(%d) ', cats{k}{1}, n);
end
fprintf('\n');

% Convert from 0-based Python to 1-based MATLAB
rxn_indices.index_matlab = rxn_indices.index_1based;

%% 7. Parameters
fprintf('[7/7] Loading all parameters...\n');
params = load_params(data_dir);

%% 8. Optional: bounds log
fpath = fullfile(data_dir, 'tracked_bounds_log.csv');
if isfile(fpath)
    bounds_log = readtable(fpath);
    fprintf('\nBounds log: %d entries\n', height(bounds_log));
else
    bounds_log = [];
end

%% Summary
fprintf('\n=== Import Complete ===\n');
fprintf('Variables in workspace:\n');
fprintf('  S            [%d x %d sparse]  - Stoichiometric matrix\n', size(S));
fprintf('  rxn_info     [%d rows]         - Reaction metadata\n', height(rxn_info));
fprintf('  met_info     [%d rows]         - Metabolite metadata\n', height(met_info));
fprintf('  lb, ub       [%d x 1]          - Flux bounds\n', length(lb));
fprintf('  rxn_indices  [%d rows]         - Categorized reaction indices\n', height(rxn_indices));
fprintf('  params       [struct]          - All kinetic/GAM/turnover parameters\n');
fprintf('  dyn_bounds   [table]           - Dynamic bound rules\n');
if ~isempty(bounds_log)
    fprintf('  bounds_log   [%d rows]         - Logged bound modifications\n', height(bounds_log));
end
fprintf('\nNext steps:\n');
fprintf('  1. Check rxn_indices for aroma/precursor reaction indices\n');
fprintf('  2. Review dyn_bounds for kinetic limit rules\n');
fprintf('  3. Use kinetic_limits_zenteno() for FBA kinetic caps\n');
fprintf('  4. Use ode_dfba_zenteno() for ODE integration\n');
fprintf('  5. Use compute_full_gam() for variable GAM\n');
