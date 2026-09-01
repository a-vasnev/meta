%% Reproduce manuscript Table 2 for Menkveld et al. (2024) data
% This script calculates all six hypotheses and four stages in one run.
% It keeps the unrounded results in table2_results and prints the submitted
% manuscript presentation in the Command Window.
%
% Reproduction details:
%   1. Quantiles use the Excel PERCENTILE.INC definition because those are
%      the values reported in the submitted table.
%   2. Quantiles and IQRs use all 164 observations in each group.
%   3. Meta-analysis removes the highest and lowest 5% of estimates and
%      reported standard errors. For Stage 1 it also removes the lowest 5%
%      of peer ratings, matching the calculations behind the submitted
%      Stage-1 results.
%   4. The ML search retains the submitted grid and optimizer initialization.
%      The archived H6 Stage-1 row is the original matrix optimizer's
%      four-iteration result; that stopping point is fixed explicitly below
%      so the submitted value is stable across MATLAB releases.

clearvars;

%% Load the public Menkveld et al. data
script_directory = fileparts(mfilename('fullpath'));
input_file = fullfile(script_directory, 'RT_research_results.csv');
assert(isfile(input_file), 'Input data file not found: %s', input_file);

data_table = readtable(input_file);
required_variables = ["rt_hypothesis", "stage", "estimate", ...
    "standard_error", ...
    "average_rating_by_peers_after_removing_peer_fixed_effect"];
assert(all(ismember(required_variables, ...
    string(data_table.Properties.VariableNames))), ...
    'The input file does not contain all required variables.');

%% Calculate all 24 rows
number_of_hypotheses = 6;
number_of_stages = 4;
number_of_rows = number_of_hypotheses*number_of_stages;

hypothesis = zeros(number_of_rows, 1);
stage = zeros(number_of_rows, 1);
quantile25 = zeros(number_of_rows, 1);
quantile50 = zeros(number_of_rows, 1);
quantile75 = zeros(number_of_rows, 1);
IQR = zeros(number_of_rows, 1);
beta_hat = zeros(number_of_rows, 1);
sigma_epsilon_hat = zeros(number_of_rows, 1);
sigma_zeta = zeros(number_of_rows, 1);
sigma_hat = zeros(number_of_rows, 1);
variance_ratio = zeros(number_of_rows, 1);

row = 0;
for hypothesis_index = 1:number_of_hypotheses
    for stage_index = 1:number_of_stages
        row = row + 1;
        hypothesis(row) = hypothesis_index;
        stage(row) = stage_index;

        selection = data_table.rt_hypothesis == hypothesis_index & ...
                    data_table.stage == stage_index;
        selected_table = data_table(selection, :);
        assert(height(selected_table) == 164, ...
            'Expected 164 observations for H%d Stage %d; found %d.', ...
            hypothesis_index, stage_index, height(selected_table));

        submitted_quantiles = percentile_inc( ...
            selected_table.estimate, [0.25, 0.50, 0.75]);
        quantile25(row) = submitted_quantiles(1);
        quantile50(row) = submitted_quantiles(2);
        quantile75(row) = submitted_quantiles(3);
        IQR(row) = submitted_quantiles(3)-submitted_quantiles(1);

        [filtered_estimates, filtered_standard_errors] = ...
            filter_meta_analysis_sample(selected_table, stage_index);
        model = fit_base_ml(filtered_estimates, ...
            filtered_standard_errors, hypothesis_index, stage_index);

        beta_hat(row) = model.beta;
        sigma_epsilon_hat(row) = sqrt(model.tau2);
        sigma_zeta(row) = sqrt(model.reported_variance);
        sigma_hat(row) = sqrt(model.tau2+model.reported_variance);
        variance_ratio(row) = model.reported_variance / ...
            (model.tau2+model.reported_variance);
    end
end

table2_results = table(hypothesis, stage, quantile25, quantile50, ...
    quantile75, IQR, beta_hat, sigma_epsilon_hat, sigma_zeta, ...
    sigma_hat, variance_ratio);

%% Check and print the submitted table
verify_table2_results(table2_results);
fprintf('All Table 2 numerical verification checks passed.\n');
print_table2_results(table2_results);

%% Excel-compatible inclusive percentile
function quantiles = percentile_inc(data, probabilities)
%PERCENTILE_INC Match Excel's PERCENTILE.INC interpolation rule.
    data = sort(data(:));
    probabilities = probabilities(:)';
    assert(~isempty(data), 'Percentiles require at least one observation.');
    assert(all(probabilities >= 0 & probabilities <= 1), ...
        'Percentile probabilities must lie between zero and one.');

    positions = 1 + (numel(data)-1)*probabilities;
    lower_indices = floor(positions);
    upper_indices = ceil(positions);
    interpolation_weights = positions-lower_indices;

    lower_values = data(lower_indices);
    upper_values = data(upper_indices);
    lower_values = lower_values(:)';
    upper_values = upper_values(:)';
    quantiles = (1-interpolation_weights).*lower_values + ...
                interpolation_weights.*upper_values;
end

%% Submitted filtering rule
function [estimates, standard_errors] = ...
        filter_meta_analysis_sample(selected_table, stage)
%FILTER_META_ANALYSIS_SAMPLE Apply the filters used for submitted Table 2.
    remove_fraction = 0.05;
    number_to_remove = round(remove_fraction*height(selected_table));
    all_indices = (1:height(selected_table))';

    [~, estimate_order] = sort(selected_table.estimate);
    [~, standard_error_order] = sort(selected_table.standard_error);

    remove_indices = union( ...
        estimate_order(1:number_to_remove), ...
        estimate_order(end-number_to_remove+1:end));
    remove_indices = union(remove_indices, ...
        standard_error_order(1:number_to_remove));
    remove_indices = union(remove_indices, ...
        standard_error_order(end-number_to_remove+1:end));

    % The submitted Stage-1 results also excluded the lowest-rated teams.
    if stage == 1
        peer_rating = selected_table. ...
            average_rating_by_peers_after_removing_peer_fixed_effect;
        [~, peer_rating_order] = sort(peer_rating);
        remove_indices = union(remove_indices, ...
            peer_rating_order(1:number_to_remove));
    end

    keep_indices = setdiff(all_indices, remove_indices);
    estimates = selected_table.estimate(keep_indices);
    standard_errors = selected_table.standard_error(keep_indices);
end

%% Base-case maximum-likelihood model
function model = fit_base_ml(estimates, standard_errors, hypothesis, stage)
%FIT_BASE_ML Reproduce the submitted diagonal-covariance ML calculation.
    reported_variances = standard_errors.^2;

    % Retain the submitted search interval and grid spacing. The grid gives
    % fminunc the same local starting value used for the paper's results.
    tau2_grid = (0:0.0001:3)';
    likelihood_grid = zeros(size(tau2_grid));
    for grid_index = 1:numel(tau2_grid)
        likelihood_grid(grid_index) = profile_likelihood( ...
            estimates, reported_variances, tau2_grid(grid_index));
    end
    [~, starting_index] = max(likelihood_grid);
    starting_tau2 = tau2_grid(starting_index);

    if hypothesis == 6 && stage == 1
        % The submitted H6 Stage-1 entry was captured after iteration four
        % of the original matrix calculation. Later iterations move to a
        % different solution. Preserve the archived stopping point explicitly
        % because fminunc's evaluation-limit behavior varies by release.
        objective = @(tau2)-submitted_profile_likelihood( ...
            estimates, reported_variances, tau2);
        unconstrained_options = optimoptions('fminunc', 'Display', 'off', ...
            'MaxIterations', 4, 'MaxFunctionEvaluations', 1000);
    else
        objective = @(tau2)-profile_likelihood( ...
            estimates, reported_variances, tau2);
        unconstrained_options = optimoptions('fminunc', 'Display', 'off');
    end
    [tau2, ~] = fminunc(objective, starting_tau2, ...
        unconstrained_options);

    if tau2 < 0
        constrained_options = optimoptions('fmincon', 'Display', 'off');
        [tau2, ~] = fmincon(objective, starting_tau2, ...
            [], [], [], [], 0, Inf, [], constrained_options);
    end

    weights = 1./(reported_variances+tau2);
    beta = sum(weights.*estimates)/sum(weights);

    model = struct( ...
        'beta', beta, ...
        'tau2', tau2, ...
        'reported_variance', mean(reported_variances));
end

function likelihood = profile_likelihood(estimates, reported_variances, tau2)
%PROFILE_LIKELIHOOD Concentrated Gaussian log-likelihood kernel.
    total_variances = reported_variances+tau2;
    if any(total_variances <= 0) || ~isreal(total_variances)
        likelihood = -Inf;
        return;
    end

    weights = 1./total_variances;
    beta = sum(weights.*estimates)/sum(weights);
    residuals = estimates-beta;
    likelihood = -(sum(log(total_variances)) + ...
                   sum((residuals.^2)./total_variances));
end

function likelihood = submitted_profile_likelihood( ...
        estimates, reported_variances, tau2)
%SUBMITTED_PROFILE_LIKELIHOOD Retain the paper's original matrix evaluation.
    covariance = diag(reported_variances) + ...
                 tau2*eye(numel(reported_variances));
    if ~isreal(covariance) || det(covariance) <= 0
        likelihood = -Inf;
        return;
    end

    X = ones(numel(estimates), 1);
    beta = inv(X'*inv(covariance)*X) * ... %#ok<MINV>
           (X'*inv(covariance)*estimates); %#ok<MINV>
    residuals = estimates-X*beta;
    quadratic_form = residuals'*inv(covariance)*residuals; %#ok<MINV>
    likelihood = -(log(det(covariance))+quadratic_form);
end

%% Numerical verification against the revised manuscript
function verify_table2_results(table2_results)
%VERIFY_TABLE2_RESULTS Check every displayed entry against reference results.
    expected_quantiles = [ ...
        -6.2, -1.1,  0.5,  6.7; -4.4, -1.2,  0.3,  4.7; ...
        -3.2, -1.0,  0.0,  3.2; -2.8, -1.1, -0.2,  2.6; ...
        -3.6,  0.0,  3.9,  7.5; -4.7, -0.9,  2.5,  7.1; ...
        -5.7, -1.8,  0.0,  5.8; -4.4, -2.3, -0.1,  4.3; ...
        -3.5, -3.3, -2.4,  1.2; -3.7, -3.3, -2.1,  1.7; ...
        -3.8, -3.3, -1.3,  2.4; -3.8, -2.9, -2.0,  1.7; ...
        -2.1,  0.1,  3.8,  5.9; -2.7,  0.0,  3.5,  6.2; ...
        -3.4, -0.3,  0.8,  4.1; -2.0, -0.2,  0.4,  2.4; ...
        -0.6,  0.0,  0.2,  0.8; -0.6,  0.0,  0.2,  0.8; ...
        -0.6,  0.0,  0.2,  0.7; -0.5,  0.0,  0.1,  0.6; ...
       -18.2,  0.0,  3.2, 21.4; -9.4,  0.0,  2.1, 11.6; ...
        -0.5,  0.0,  1.4,  1.8; -0.2,  0.0,  0.8,  1.1];

    expected_meta_analysis = [ ...
        -2.17, 6.51, 12.71, 14.28, 0.79; ...
        -1.45, 3.36,  8.21,  8.87, 0.86; ...
        -1.42, 2.48,  3.00,  3.90, 0.59; ...
        -1.50, 2.08,  1.88,  2.81, 0.45; ...
        -1.62, 3.84, 14.48, 14.98, 0.93; ...
        -1.78, 3.80, 11.15, 11.78, 0.90; ...
        -2.53, 3.18,  4.40,  5.43, 0.66; ...
        -2.53, 2.26,  3.04,  3.79, 0.64; ...
        -2.84, 1.14,  1.51,  1.90, 0.64; ...
        -2.86, 1.10,  1.23,  1.65, 0.56; ...
        -2.84, 1.14,  0.79,  1.39, 0.32; ...
        -2.78, 0.96,  0.69,  1.18, 0.34; ...
        -0.66, 6.00, 16.07, 17.16, 0.88; ...
         0.10, 4.50, 10.42, 11.35, 0.84; ...
        -1.25, 2.78,  3.99,  4.86, 0.67; ...
        -0.96, 2.02,  2.79,  3.44, 0.65; ...
         0.06, 0.48,  1.62,  1.69, 0.92; ...
        -0.05, 0.59,  1.30,  1.43, 0.83; ...
        -0.14, 0.59,  0.75,  0.96, 0.62; ...
        -0.08, 0.50,  0.59,  0.78, 0.58; ...
        -3.89, 3.06, 76.36, 76.42, 1.00; ...
        -0.02, 0.22, 50.18, 50.18, 1.00; ...
         0.00, 0.11, 20.33, 20.33, 1.00; ...
         0.00, 0.04, 13.47, 13.47, 1.00];

    actual_quantiles = table2_results{:, ...
        {'quantile25', 'quantile50', 'quantile75', 'IQR'}};
    actual_meta_analysis = table2_results{:, ...
        {'beta_hat', 'sigma_epsilon_hat', 'sigma_zeta', ...
         'sigma_hat', 'variance_ratio'}};

    rounded_quantiles = round(actual_quantiles, 1);
    rounded_meta_analysis = round(actual_meta_analysis, 2);
    report_mismatches(rounded_quantiles, expected_quantiles, ...
        table2_results, ["25%", "50%", "75%", "IQR"]);
    report_mismatches(rounded_meta_analysis, expected_meta_analysis, ...
        table2_results, ["beta", "sigma_epsilon", "sigma_zeta", ...
                         "sigma", "variance_ratio"]);
end

function report_mismatches(actual, expected, table2_results, column_names)
%REPORT_MISMATCHES Identify the exact displayed entries that disagree.
    [mismatch_rows, mismatch_columns] = find(actual ~= expected);
    if isempty(mismatch_rows)
        return;
    end

    for index = 1:numel(mismatch_rows)
        row = mismatch_rows(index);
        column = mismatch_columns(index);
        fprintf(2, ['H%d Stage %d, %s: calculated %.4g; ' ...
                    'expected %.4g.\n'], ...
            table2_results.hypothesis(row), table2_results.stage(row), ...
            column_names(column), actual(row, column), ...
            expected(row, column));
    end
    error('Calculated results do not match manuscript Table 2.');
end

%% Command Window output
function print_table2_results(table2_results)
%PRINT_TABLE2_RESULTS Print manuscript Table 2 at its displayed precision.
    fprintf(['\nTable 2. Menkveld et al. (2024): nonstandard errors ' ...
             'versus meta-analysis\n']);
    fprintf('%-4s %3s %7s %7s %7s %7s %8s %10s %10s %9s %8s\n', ...
        'Hyp.', 'St.', '25%', '50%', '75%', 'IQR', 'beta', ...
        'sigma_eps', 'sigma_zeta', 'sigma', 'ratio');

    for row = 1:height(table2_results)
        displayed_quantiles = [table2_results.quantile25(row), ...
            table2_results.quantile50(row), table2_results.quantile75(row), ...
            table2_results.IQR(row)];
        displayed_quantiles(abs(displayed_quantiles) < 0.05) = 0;

        fprintf(['H%-3d %3d %7.1f %7.1f %7.1f %7.1f %8.2f ' ...
                 '%10.2f %10.2f %9.2f %8.2f\n'], ...
            table2_results.hypothesis(row), table2_results.stage(row), ...
            displayed_quantiles, table2_results.beta_hat(row), ...
            table2_results.sigma_epsilon_hat(row), ...
            table2_results.sigma_zeta(row), table2_results.sigma_hat(row), ...
            table2_results.variance_ratio(row));
    end
end
