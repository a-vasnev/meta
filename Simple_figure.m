% Create a new figure
close all; clear all;
%figure;
scrsz = get(0,'ScreenSize'); gcf=figure('PaperPositionMode','auto','Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/(4.5*1.1) scrsz(4)/(4.5*1.1)],'PaperOrientation','landscape'); hold on; % position [left bottom width height]

x1 = 72;
sigma1 = 1;
sigma2 = 1;

% Set x-axis limits
%xlim([x1-3, x1+3]);
xlim([x1-5, x1+5]);

% Set y-axis limits (change this to fit your needs)
ylim([0, 2.5]);

% Draw the horizontal line at y = 1/2
hold on; % Allows multiple plots on the same figure
line([x1-2*sigma1, x1+2*sigma1], (sigma1^2 + sigma2^2)/4*[1, 1], 'Color', 'b', 'LineStyle', '--');

% Draw the horizontal line at y = 1
line([0, 100], sigma1^2 * [1, 1], 'Color', 'r', 'LineStyle', '--');

% Draw the vertical line at x1
line(x1*[1, 1], [0, 2.5], 'Color', 'r', 'LineStyle', ':');

% Constants for the hyperbolic function
b1 = 0.5*(72-0)*(72-72*2);
b2 = 0.5*(72-44)*(72-100);

% Generate x-values, avoiding the exact points where the function would be undefined
x_values = linspace(x1-10, x1+10, 1000)';

% Compute y-values for the hyperbolic curve
%y_values = (b1 ./ ((x_values - 0) .* (x_values - 72*2))) .* (x_values <= 72) + (b2 ./ ((x_values - 44) .* (x_values - 100))) .* (x_values > 72);
%y_values =  (b2 ./ ((x_values - 44) .* (x_values - 100))) ;

% do parabola instead
tmp = (1-0.5) / ((x1-2*sigma1-x1)^2);
y_values = tmp * (x_values - x1) .^ 2 + 0.5;

% ML curve
y_values_ML = 1/8 * (x_values - x1) .^ 2; 
%y_values_ML = ((x_values-x1) >= 2) *  1/8 * (x_values - x1) .^ 2 + ((x_values-x1) < 2) * 0;
for j = 1:length(x_values)
   if abs(x_values(j)-x1) < 2
       %y_values_ML(j) = nan;
   end
end

% stat framework
%y_values_stat = cdf('Normal',x1,x_values,sigma1) ./ (1 - cdf('Normal',x1,x_values,sigma1));
y_values_stat = (1-cdf('Normal',x1,x_values,sigma1)) .* (x_values <= x1) + (cdf('Normal',x1,x_values,sigma1)) .* (x_values > x1);
%y_values_stat = 0.01*(1./(1-y_values_stat) - 2) + 0.5; % transform that it goes to infinity

% meta-analysis random-effects model

y_values_meta = y_values_stat * 0;
k = 2; % number of studies
df = k -1;
w1 = 1/(sigma1^2);
w2 = 1/(sigma2^2);
C = w1 + w2 - (w1^2 + w2^2)/(w1 + w2);
for j = 1:length(x_values)
    x2 = x_values(j);
    Q = w1 * x1^2 + w2 * x2^2 - (w1 * x1 + w2 * x2)^2 / (w1 + w2);
    T_sq = (Q - df) / C;
    if T_sq > 0
        w1_star = 1 / (sigma1^2 + T_sq);
        w2_star = 1 / (sigma2^2 + T_sq);
        M_star = (w1_star * x1 + w2_star * x2) / (w1_star + w2_star);
        VM_star = 1 / (w1_star + w2_star);
        y_values_meta(j) = VM_star;
    else
        y_values_meta(j) = nan;
    end
end

% Plot the curve
plot(x_values, y_values, 'g');
%plot(x_values, y_values_stat, 'r');
plot(x_values, y_values_meta, 'k');
plot(x_values, y_values_ML, 'b');
plot(72,1/2,'or')

% Add labels and title
xlabel('y_2');
ylabel('variance of combination');
%title('Simple example: x_1 = 72, \sigma_1=\sigma_2=1');

% Add legend
%legend('y = 1/2 (stat)', 'y = 1 (ignore x_2)', 'x_2=x_1', 'Parabola','stat compatability','meta-analysis RE','min','Location','southeast');
%legend('y = 1/2 (stat)', 'y = 1 (ignore x_2)', 'x_2=x_1', 'Parabola','meta-analysis RE','min','Location','southeast');
%legend('', '', '', 'Parabola','Method of moments','Maximum likelihood','Location','northeast');
legend('', '', '', 'Common sense','MOM','ML','Location','northeast');

hold off; % Stops adding to the current figure

% print into file
print_flag = 0;
if (print_flag == 1) 
    path_name = ['/Users/av/Dropbox (Sydney Uni)/Research/2023/IPCC/Draft_meta/']; % mac path name
    %path_name = ['C:\Users\avasnev\Dropbox (Sydney Uni)\Research\2023\IPCC\Draft_meta\']; % PC path name
    %path_name = []; % online
    output_name = [path_name 'Figure_meta_simple_figure.eps'];
    print(output_name,'-dpsc2');
    %print(output_name,'-depsc'); % formats '-dpsc2' - I used before; '-dpdf' - doesn't have bounding box
    %exportgraphics(gcf,output_name,'BackgroundColor','none','ContentType','vector')
end 
