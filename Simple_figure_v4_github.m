%% Simple two-input figures
close all;
clearvars;

%% Configuration
% Set each element to true to export the corresponding manuscript figure
% as a vector PDF: [Figure 1, Figure 2].
export_figures = [false, false];

script_directory = fileparts(mfilename('fullpath'));
output_directory = fullfile(script_directory, 'output', 'pdf');
if any(export_figures) && ~isfolder(output_directory)
    mkdir(output_directory);
end

%% Inputs
x1 = 72;
sigma1 = 1;
sigma2 = 1;
minimum_variance = (sigma1^2 + sigma2^2)/4;

%% Variance curves
x_values = linspace(x1-10, x1+10, 1000)';

% Common-sense parabola: its minimum is 1/2 at y_2 = y_1, and it reaches
% one at y_2 = y_1 +/- 2*sigma1.
common_sense_scale = (1-0.5)/((x1-2*sigma1-x1)^2);
variance_common_sense = common_sense_scale*(x_values-x1).^2 + 0.5;

% Maximum-likelihood curve. Within two standard deviations of x1, its
% variance is the constant 1/2 shown by the dashed blue line.
variance_ml = (1/8)*(x_values-x1).^2;
variance_ml(abs(x_values-x1) < 2) = nan;

% Method-of-moments random-effects variance.
variance_mom = nan(size(x_values));
number_of_studies = 2;
degrees_of_freedom = number_of_studies - 1;
weight1 = 1/sigma1^2;
weight2 = 1/sigma2^2;
weight_adjustment = weight1 + weight2 ...
                  - (weight1^2 + weight2^2)/(weight1 + weight2);

for j = 1:length(x_values)
    x2 = x_values(j);
    q_statistic = weight1*x1^2 + weight2*x2^2 ...
                - (weight1*x1 + weight2*x2)^2/(weight1 + weight2);
    between_study_variance = ...
        (q_statistic-degrees_of_freedom)/weight_adjustment;

    if between_study_variance > 0
        random_weight1 = 1/(sigma1^2 + between_study_variance);
        random_weight2 = 1/(sigma2^2 + between_study_variance);
        variance_mom(j) = 1/(random_weight1 + random_weight2);
    end
end

%% Figure 1: common-sense parabola
% Fixed dimensions keep the layout consistent across computers.
common_sense_figure = figure( ...
    'Units', 'pixels', ...
    'PaperPositionMode', 'auto', ...
    'Position', [100, 100, 540, 372], ...
    'PaperOrientation', 'landscape', ...
    'Name', 'Common-sense variance', ...
    'NumberTitle', 'off');
hold on;

xlim([x1-4, x1+4]);
ylim([0, 2.5]);
draw_reference_lines(x1, sigma1, minimum_variance);
common_sense_curve = plot(x_values, variance_common_sense, 'g', ...
                          'LineWidth', 1.5);
plot(x1, minimum_variance, 'or', 'HandleVisibility', 'off');

xlabel('y_2');
ylabel('variance of combination');
legend(common_sense_curve, {'Common sense'}, 'Location', 'northeast');
hold off;

%% Figure 2: common-sense, MOM, and ML variances
combined_figure = figure( ...
    'Units', 'pixels', ...
    'PaperPositionMode', 'auto', ...
    'Position', [100, 100, 550, 372], ...
    'PaperOrientation', 'landscape', ...
    'Name', 'Random-effects meta-analysis', ...
    'NumberTitle', 'off');
hold on;

xlim([x1-3, x1+3]);
ylim([0, 2.5]);
draw_reference_lines(x1, sigma1, minimum_variance);

common_sense_curve = plot(x_values, variance_common_sense, 'g');
mom_curve = plot(x_values, variance_mom, 'k');
ml_curve = plot(x_values, variance_ml, 'b');
plot(x1, minimum_variance, 'or', 'HandleVisibility', 'off');

xlabel('y_2');
ylabel('variance of combination');
legend([common_sense_curve, mom_curve, ml_curve], ...
       {'Common sense', 'MOM', 'ML'}, ...
       'Location', 'northeast');
hold off;

%% Optional PDF export
if export_figures(1)
    exportgraphics(common_sense_figure, ...
        fullfile(output_directory, ...
                 'Figure_meta_simple_figure_common_sense.pdf'), ...
        'BackgroundColor', 'none', ...
        'ContentType', 'vector');
end

if export_figures(2)
    exportgraphics(combined_figure, ...
        fullfile(output_directory, 'Figure_meta_simple_figure.pdf'), ...
        'BackgroundColor', 'none', ...
        'ContentType', 'vector');
end

%% Draw reference lines shared by both figures
function draw_reference_lines(x1, sigma1, minimum_variance)
%DRAW_REFERENCE_LINES Add the equal-weight, single-input, and agreement lines.
    line([x1-2*sigma1, x1+2*sigma1], ...
         minimum_variance*[1, 1], ...
         'Color', 'b', ...
         'LineStyle', '--', ...
         'HandleVisibility', 'off');

    line([0, 100], sigma1^2*[1, 1], ...
         'Color', 'r', ...
         'LineStyle', '--', ...
         'HandleVisibility', 'off');

    line(x1*[1, 1], [0, 2.5], ...
         'Color', 'r', ...
         'LineStyle', ':', ...
         'HandleVisibility', 'off');
end
