%% v5 added Menkveld et al (2024) data as my_case = 12
%% v6_github keeps only cases used in the paper

%% clean up
close all; clear all;


%% define constants and data
my_q = norminv(0.975); % normal quantile to replace 1.96

my_case = 12; % selects the appropriate data below; use 1,2,4,5,6 for now, the rest needs more work
             % 2 - corresponds to Doi-1 case in our paper
             % 12 - Menkveld et al (2024) data

if my_case == 2
    % data from Doi et al (2015-I) fig 3 based on Collins et al (1985) - not symmetric and bounded by 0
    Y_i =  [1.04 0.40 0.33 0.23 0.25 0.74 0.77 2.97 1.14]'; % point estimates
    LB_i = [0.48 0.20 0.14 0.08 0.13 0.59 0.39 0.59 0.69]'; % lower bound of 95% conf interval
    UB_i = [2.28 0.78 0.74 0.67 0.48 0.94 1.52 15.07 1.91]'; % upper bound of 95% conf interval
    % it's symmetric in logs, but rounding causes discreptancies, the
    % largest difference between left and right is 0.05 for observation 3
    % the original paper doesn't have the ratios as in Doi et all
    % decision: use Y_i and LB_i
    Y_i = log(Y_i);
    LB_i = log(LB_i);
    UB_i = log(UB_i);
    my_title = 'Example from Doi et al (2015-I) fig 3 based on Collins et al (1985) - IN LOGS';
    my_dates = [1962 1962 1964 1964 1964 1965 1966 1971 1975]';


elseif my_case == 12 % Menkveld et al (2024) data; choose Hypothesis and Stage here
    % Define the full path
    filePath = '/Users/avasnev/Sydney Uni Dropbox/Andrey Vasnev/Research/2023/IPCC/Figures/Nonstandard errors/AV analysis/RT_research_results.csv';
    % Read the CSV file into a table
    dataTable = readtable(filePath);
    % get quantiles accross all RT-H and all stages
    summaryStats_all = [];
    for my_hypothesis = 1:6
        for my_stage = 1:4
            my_selection = (dataTable.rt_hypothesis == my_hypothesis)&(dataTable.stage == my_stage);
            dataTable_filtered = dataTable(my_selection,:);
            Y_i = dataTable_filtered.estimate;
            summaryStats_all = [summaryStats_all;
                quantile(Y_i, [0 0.10 0.25 0.50 0.75 0.90 1], Method="approximate")];
        end
    end
    % zoom in on particular RT-H and stage
    my_hypothesis = 1;
    my_stage = 3; % focus only on stage 1 for now
    % select slice of the data for one hypothesis and one stage
    my_selection = (dataTable.rt_hypothesis == my_hypothesis)&(dataTable.stage == my_stage);
    dataTable_filtered = dataTable(my_selection,:);
    if my_hypothesis == 0 % switch this case off, so the same filter is applied to all hypotheses
        my_points = [1:27 29:164]; % remove CXHI outlier for RT-H1, stage 1
        %my_points = [1:164]; % all points
    elseif my_hypothesis >= 1
        %my_points = [1:108 110:143 145:164]; % remove PVCC, VHSY outliers for RT-H2, stage 1
        % or alternatively remove top and bottom 5%
        Y_i = dataTable_filtered.estimate;
        remove_pct = 0.05; % percentage to remove
        n_extr = round(remove_pct * length(Y_i)); % number of extreme observations to remove
        % Sort and get original indices
        [~, sorted_idx] = sort(Y_i);
        % Indices of 5% smallest
        idx_smallest_5pct = sorted_idx(1:n_extr);
        % Indices of 5% largest
        idx_largest_5pct = sorted_idx(end-n_extr+1:end);
        % also remove outliers in st deviation
        sigma_i = dataTable_filtered.standard_error;
        [~, sorted_idx_sigma] = sort(sigma_i);
        idx_smallest_5pct_sigma = sorted_idx_sigma(1:n_extr);
        idx_largest_5pct_sigma = sorted_idx_sigma(end-n_extr+1:end);
        all_idx = [1:length(Y_i)]';
        remove_idx = union(idx_smallest_5pct, idx_largest_5pct);  % combine both sets
        remove_idx = union(remove_idx, idx_largest_5pct_sigma);  % add outliers in st dev
        remove_idx = union(remove_idx, idx_smallest_5pct_sigma);
        % add team quality to the filter
        if my_stage ==1
            peer_rating_i = dataTable_filtered.average_rating_by_peers_after_removing_peer_fixed_effect;
            [~, sorted_idx_peer_rating] = sort(peer_rating_i);
            idx_smallest_5pct_peer_rating = sorted_idx_peer_rating(1:n_extr);
            remove_idx = union(remove_idx, idx_smallest_5pct_peer_rating);
        end
        % alternatively create Loss of Confidence Set
        %remove_idx = remove_LCS(dataTable); % it removes 125 out of 164 teams!!!
        my_points = setdiff(all_idx, remove_idx);  % keep the rest
    else
        my_points = [1:164];
    end
    %figure; plot(Y_i, 1./sigma_i,'.'); hold on; % funnel plot on unfiltered data
    Y_i = dataTable_filtered.estimate(my_points);
    sigma_i = dataTable_filtered.standard_error(my_points);
    my_title = ['Menkveld et al (2024) data; RT-H', num2str(my_hypothesis), '; stage ', num2str(my_stage)];
    %figure; plot(Y_i, 1./sigma_i,'.'); hold on; % funnel plot on FILTERED data
    % kernel density plot for filtered data
    %[f1, xi1] = ksdensity(Y_i);
    %figure; plot(xi1, f1, 'LineWidth', 2); hold on; title("Kernel density of Y_i")
else
    error("no such case yet")
end

if (my_case <= 4) || ((my_case >= 9) && (my_case <= 11))
    sigma_i = (Y_i - LB_i)/my_q; % uncertainty of individual studies
else
    LB_i = Y_i - my_q * sigma_i; % Lower bounds for IPCC studies and Mankveld et al (2024) data
end

sigma2_i = sigma_i.^2;
w_i = 1./(sigma2_i);

%% fixed effect as in Borenstein et al (2009)
M_FE = w_i'*Y_i / sum(w_i);
V_FE = 1 / sum(w_i);

%% random effect as in Borenstein et al (2009)
k = length(Y_i); % number of studies
df = k -1;
C = sum(w_i) - w_i'*w_i / sum(w_i);
Q = w_i'*(Y_i.^2) - ((w_i'*Y_i)^2) / sum(w_i); 
T2 = (Q - df) / C;

w_i_star = 1./(sigma2_i + T2);

M_RE = w_i_star'*Y_i / sum(w_i_star);
V_RE = 1 / sum(w_i_star);

CI_RE = [M_RE - my_q * sqrt(V_RE); M_RE + my_q * sqrt(V_RE)]; 

%% OLS as per Jan's note from mata-analysis-reflections.pdf

M_OLS = mean(Y_i);
sigma2_eps_OLS = mean((Y_i - M_OLS).^2);
CI_OLS = [M_OLS - my_q * sqrt(sigma2_eps_OLS/k); M_OLS + my_q * sqrt(sigma2_eps_OLS/k)];

%% Our analysis

V = diag(sigma2_i);

X = ones(k,1);

tau2_set = [0.0:0.0001:3]';
phi_set =  tau2_set * 0;
for j = 1:length(tau2_set)
   phi_set(j) = phi_tau(Y_i, X, V, tau2_set(j));
end
%plot(tau2_set, phi_set) % graph: concentrated likelihood as function of tau

[tmp,I] = max(phi_set);
tau2_ML_grid = tau2_set(I);
fun = @(tau2)-phi_tau(Y_i, X, V, tau2);
[tau2_ML, fval] = fminunc(fun,tau2_ML_grid);
if tau2_ML < 0 % add non negativity constraint
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [0];
    ub = [Inf];
    [tau2_ML, fval] = fmincon(fun,tau2_ML_grid,A,b,Aeq,beq,lb,ub);
end

M_ML = inv(X' * inv(Omega(V,tau2_ML)) * X) * (X' * inv(Omega(V,tau2_ML)) * Y_i);

%% figure; 

if my_case == 8
    x_text = min(LB_i) - 5.1;
    y_text_delta = 0.25;
    scrsz = get(0,'ScreenSize'); gcf=figure('PaperPositionMode','auto','Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/(4.5) scrsz(4)/(6)],'PaperOrientation','landscape'); hold on; % position [left bottom width height]
elseif my_case == 2
    x_text = min(LB_i)-0.3;
    y_text_delta = 0.4;
    scrsz = get(0,'ScreenSize'); gcf=figure('PaperPositionMode','auto','Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/(4.5) scrsz(4)/(2.5)],'PaperOrientation','landscape'); hold on; % position [left bottom width height]
elseif (my_case >= 9) && (my_case <=11)
    x_text = min(LB_i)-0.3;
    y_text_delta = 0.4;
    scrsz = get(0,'ScreenSize'); gcf=figure('PaperPositionMode','auto','Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/(4.5) scrsz(4)/(2.5)],'PaperOrientation','landscape'); hold on; % position [left bottom width height]
elseif (my_case == 1) || (my_case == 6)
    x_text = min(LB_i) + 0.007*(my_case == 1) + 0.05*(my_case == 6);
    y_text_delta = 0;
    scrsz = get(0,'ScreenSize'); gcf=figure('PaperPositionMode','auto','Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/(4.5) scrsz(4)/(3.5)],'PaperOrientation','landscape'); hold on; % position [left bottom width height]
elseif (my_case == 4) || (my_case == 5)
    x_text = min(LB_i) - 0.45*(my_case == 4) - 0.15*(my_case == 5);
    y_text_delta = 0;
    scrsz = get(0,'ScreenSize'); gcf=figure('PaperPositionMode','auto','Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/(4.5) scrsz(4)/(2.5)],'PaperOrientation','landscape'); hold on; % position [left bottom width height]
elseif (my_case == 12)
    x_text = min(LB_i)-0.3;
    y_text_delta = 0.4;
    scrsz = get(0,'ScreenSize'); 
    % to view on screen (my MacBook)
    %gcf=figure('PaperPositionMode','auto','Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/(2.0) scrsz(4)/(1.1)],'PaperOrientation','landscape'); hold on; % position [left bottom width height]
    % to print for latex file
    gcf=figure('PaperPositionMode','auto','Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/(2.0*1.5) scrsz(4)/(1.1*1.5)],'PaperOrientation','portrait'); hold on; % position [left bottom width height]
    box on;
else    
    x_text = min(LB_i);
    y_text_delta = 0;
    scrsz = get(0,'ScreenSize'); gcf=figure('PaperPositionMode','auto','Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/(4.5) scrsz(4)/(3.5)],'PaperOrientation','landscape'); hold on; % position [left bottom width height]
end

if (my_case < 9)
    plot(Y_i,[k:-1:1]','bs'); hold on; % individual studies
    for i = 1:k
        plot([Y_i(i) - my_q*sigma_i(i); Y_i(i) + my_q*sigma_i(i)],(k-i+1)*[1;1],'b-')
        text(x_text, k-i+1,num2str(i))
    end
elseif ((my_case >= 9) && (my_case <= 11))
    k_case2 = 9;
    for i = 1:k
        plot(Y_i(i),(k_case2-my_points(i)+1)','bs'); hold on; % individual studies
        plot([Y_i(i) - my_q*sigma_i(i); Y_i(i) + my_q*sigma_i(i)],(k_case2-my_points(i)+1)*[1;1],'b-')
        text(x_text, k_case2-my_points(i)+1,num2str(my_points(i)))
    end
elseif (my_case == 12)
    k_case2 = 164;
    for i = 1:k
        plot(Y_i(i),(k_case2-my_points(i)+1)','bs'); hold on; % individual studies
        plot([Y_i(i) - my_q*sigma_i(i); Y_i(i) + my_q*sigma_i(i)],(k_case2-my_points(i)+1)*[1;1],'b-')
        %text(x_text, k_case2-my_points(i)+1,num2str(my_points(i)))
    end
else
    errror("need to adjust figure cases")
end

%plot(M_FE, -1,'rd'); % fixed effect
%plot([M_FE - my_q*sqrt(V_FE); M_FE + my_q*sqrt(V_FE)],(-1)*[1;1],'r-')
%%plot([M_FE - my_q*sqrt(V_FE); M_FE + my_q*sqrt(V_FE)],(-1)*[1;1],'k|')
%text(x_text, -1,'FE')

%plot([M_RE - my_q*sqrt(V_RE); M_RE + my_q*sqrt(V_RE)],(-2)*[1;1],'k|')
if (my_case == 1) || (my_case == 2) || (my_case == 4) || (my_case >= 9)
    %plot(M_RE, -2,'rd'); % random effect
    %plot([M_RE - my_q*sqrt(V_RE); M_RE + my_q*sqrt(V_RE)],(-2)*[1;1],'r-'); 
    %text(x_text, -2+y_text_delta/2,'RE')
elseif my_case == 8
    plot(M_RE, -3,'rd'); % random effect
    plot([M_RE - my_q*sqrt(V_RE); M_RE + my_q*sqrt(V_RE)],(-3)*[1;1],'r-'); 
    text(x_text, -3+y_text_delta,'prediction')
else
    plot(M_RE, -2,'rd'); % random effect
    plot([M_RE - my_q*sqrt(V_RE); M_RE + my_q*sqrt(V_RE)],(-2)*[1;1],'r-'); 
    text(x_text, -2+y_text_delta,'RE')
end

if my_case ~= 8
    if (my_case == 1) || (my_case == 2) || (my_case == 4) || (my_case >= 9)
        position = -1; % y-coordinate for the line and text
    else
        position = -3;
    end
    plot(M_ML, position,'rd'); % ML - base case
    V_ML_y = tau2_ML + trace(V)/k;
    st_err_ML = sqrt(V_ML_y);
    V_ML = inv(X' * inv(Omega(V,tau2_ML)) * X);
    st_err_M_ML = sqrt(V_ML);
    plot([M_ML - my_q*st_err_ML; M_ML + my_q*st_err_ML],(position)*[1;1],'r-') % full interval with tau2_ML + trace(V)/k
    %plot([M_ML - my_q*sqrt(trace(V)/k); M_ML + my_q*sqrt(trace(V)/k)],(-3)*[1;1],'kx') % interval with only trace(V)/k
    plot([M_ML - my_q*st_err_M_ML; M_ML + my_q*st_err_M_ML],(position)*[1;1],'k|') % interval for mu
    text(x_text, position+y_text_delta,'Base case')
    my_output_table4 = [M_ML sqrt(tau2_ML) sqrt(trace(V)/k) st_err_ML (trace(V)/k)/V_ML_y];
    my_output_fig7 = [M_ML, M_ML - my_q*st_err_ML, M_ML + my_q*st_err_ML, M_ML - my_q*st_err_M_ML, M_ML + my_q*st_err_M_ML];
    my_q_50p = norminv(0.75);
    my_output_fig7_50p = [M_ML, M_ML - my_q_50p*st_err_ML, M_ML + my_q_50p*st_err_ML, M_ML - my_q_50p*st_err_M_ML, M_ML + my_q_50p*st_err_M_ML];
elseif my_case == 8
    plot(M_OLS, -2,'rd'); % OLS
    plot(CI_OLS,(-2)*[1;1],'r-')
    text(x_text, -2,'RE')
end


if my_case == 7 
    plot(mu_zero, 0,'gd'); % true parameters
    st_err_zero = sqrt(tau2_zero + trace(V)/k);
    st_err_M_zero = sqrt(inv(X' * inv(Omega(V,tau2_zero)) * X));
    plot([mu_zero - my_q*st_err_zero; mu_zero + my_q*st_err_zero],(0)*[1;1],'g-') % full interval with tau2_ML + trace(V)/k
    plot([mu_zero - my_q*sqrt(trace(V)/k); mu_zero + my_q*sqrt(trace(V)/k)],(0)*[1;1],'kx') % interval with only trace(V)/k
    plot([mu_zero - my_q*st_err_M_zero; mu_zero + my_q*st_err_M_zero],(0)*[1;1],'k|') % interval for mu
    text(x_text, 0,'True')
end

if my_case < 9
    plot(M_ML*[1 1], [-6 k+0.5], ':k') % vertical line at M_ML
else
    plot(M_ML*[1 1], [-6 k_case2+0.5], ':k') % vertical line at M_ML
end

if my_case == 8
    ylim([-3.5 k+0.5])
elseif (my_case == 6) 
    ylim([-4.5 k+0.5])
elseif (my_case == 1) || (my_case == 2) || (my_case == 4)
    %ylim([-3.5 k+0.5])
    ylim([-4 k+0.5])
elseif (my_case >= 9) && (my_case <= 11)
    ylim([-4 k_case2+0.5])
    xlim([-3 3])
elseif (my_case == 12)
    ylim([-10 k_case2 + 2])
    %xlim([-600 600])
    %title(my_title) % remove title when producing files for latex
else
    ylim([-5.5 k+0.5])
end
%title(my_title)
yticklabels({}) % remove default axis labels

%% ML - case 3: relative model
fun3 = @(x)-phi_tau(Y_i, X, exp(x(2))*V, exp(x(1))); % tau2 = exp(x(1)), sigma2_eps = exp(x(2))
if my_case<6
    x0 = [log(tau2_ML), 0];
else
    x0 = [log(0.01), 0];
end
[x_ML3, fval] = fminunc(fun3,x0);
tau2_ML3 = exp(x_ML3(1));
sigma2_eps = exp(x_ML3(2));
M_ML3 = inv(X' * inv(Omega(sigma2_eps*V,tau2_ML3)) * X) * (X' * inv(Omega(sigma2_eps*V,tau2_ML3)) * Y_i);
V_ML3_y = tau2_ML3 + trace(sigma2_eps*V)/k;
st_err_ML3 = sqrt(V_ML3_y);
V_ML3 = inv(X' * inv(Omega(sigma2_eps*V,tau2_ML3)) * X);
st_err_M_ML3 = sqrt(V_ML3);
if my_case ~= 8
    if (my_case == 1) || (my_case == 2) || (my_case == 4) || ((my_case >= 9) && (my_case <= 11))
        position = -2;
    elseif (my_case == 12)
        position = -5;
    else
        position = -4;
    end
    if my_case ~= 12
        plot(M_ML3, position,'rd'); % ML - case 3
        plot([M_ML3 - my_q*st_err_ML3; M_ML3 + my_q*st_err_ML3],(position)*[1;1],'r-') % full interval with tau2_ML + trace(V)/k
        %plot([M_ML3 - my_q*sqrt(trace(sigma2_eps*V)/k); M_ML3 + my_q*sqrt(trace(sigma2_eps*V)/k)],(-4)*[1;1],'kx') % interval with only trace(V)/k
        plot([M_ML3 - my_q*st_err_M_ML3; M_ML3 + my_q*st_err_M_ML3],(position)*[1;1],'k|') % interval for mu
        text(x_text, position+y_text_delta,'Relative case')
    end
end
% % look at the likelihood surfice 
% figure
% [tau2_mesh,sigma2_mesh] = meshgrid(0.05:.05:5);
% tau2_mesh = tau2_ML3 * tau2_mesh;
% sigma2_mesh = 1 * sigma2_mesh;
% phi_ML3_mesh = tau2_mesh * 0; % preallocation
% for i = 1:size(tau2_mesh,1)
%     for j = 1:size(tau2_mesh,2)
%         phi_ML3_mesh(i,j) = phi_tau(Y_i, X, sigma2_mesh(i,j)*V, tau2_mesh(i,j));
%     end
% end
% mesh(tau2_mesh,sigma2_mesh,phi_ML3_mesh)
% xlabel('\tau^2')
% ylabel('\sigma^2')
% %xlim([0 0.0025]); ylim([0 2.5]);  
% if my_case == 1
%     zlim([30 45]);
% elseif my_case == 2
%     zlim([-5 0]);
% elseif my_case == 4
%     zlim([30 45]);
% elseif my_case == 5
%     zlim([10 40]);
% end
% % profile likelihood as function of tau^2
% figure
% subplot(1,2,1)
% tau2_set = [0.0001:0.0001:0.035]';
% phi_set =  tau2_set * 0;
% for j = 1:length(tau2_set)
%    phi_set(j) = phi_tau(Y_i, X, sigma2_eps*V, tau2_set(j));
% end
% plot(tau2_set, phi_set) % graph: concentrated likelihood as function of tau2
% hold on
% plot(tau2_ML3,phi_tau(Y_i, X, sigma2_eps*V, tau2_ML3),'ro')
% if my_case == 1
%     ylim([42 45])
% elseif my_case == 2
%     ylim([-1 0])
% elseif my_case == 4
%     ylim([40 41])
% elseif my_case == 5
%     ylim([30 40])
% elseif my_case == 6
%     ylim([0 20])
% end
% xlabel('\tau^2')
% % profile likelihood as function of sigma^2
% subplot(1,2,2)
% sigma2_set = [0:0.0001:5]';
% phi_set =  sigma2_set * 0;
% for j = 1:length(sigma2_set)
%    phi_set(j) = phi_tau(Y_i, X, sigma2_set(j)*V, tau2_ML3);
% end
% plot(sigma2_set, phi_set) % graph: concentrated likelihood as function of tau2
% hold on
% plot(sigma2_eps,phi_tau(Y_i, X, sigma2_eps*V, tau2_ML3),'ro')
% if my_case == 1
%     ylim([42 45])
% elseif my_case == 2
%     ylim([-1 0])
% elseif my_case == 4
%     ylim([40 41])
% elseif my_case == 5
%     ylim([30 40])
% elseif my_case == 6
%     ylim([0 20])
% end
% xlabel('\sigma^2')

if my_case == 12
    %summaryStats = prctile(Y_i, [25 50 75]);
    Y_i_original = dataTable_filtered.estimate;
    summaryStats = prctile(Y_i_original, [25 50 75]); % use quantiles from the full sample
    position = -8;
    q3 = norminv(0.75);
    plot(summaryStats(2), position,'rd'); % Median
    % stretch NSE interval
    %plot([summaryStats(1); summaryStats(3)],(position)*[1;1],'kx'); % IQR
    plot([summaryStats(2)+(summaryStats(1)-summaryStats(2))*my_q/q3;...
        summaryStats(2)+(summaryStats(3)-summaryStats(2))*my_q/q3],(position)*[1;1],'r-') % NSE
    text(x_text, position+y_text_delta,'NSE interval (adjusted)')
    % check stats with outlier; they don't match Table 1 in the paper!!!
    % They match with quantile in Excel so the difference is in quantile
    % function
    
    % figure with two intervals
    my_output_fig7 = [my_output_fig7, summaryStats];
    my_output_fig7_50p = [my_output_fig7_50p, summaryStats];
    % stretch NSE interval
    q3 = norminv(0.75);
    my_output_fig7(5+1) = my_output_fig7(5+2)+(my_output_fig7(5+1) - my_output_fig7(5+2))*my_q/q3;
    my_output_fig7(5+3) = my_output_fig7(5+2)+(my_output_fig7(5+3) - my_output_fig7(5+2))*my_q/q3;
    x_text = my_output_fig7(2);
    if (1==0) % switch it off to print the big picture into file
        gcf2 = figure('PaperPositionMode','auto','Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/(5.0) scrsz(4)/(7.5)],'PaperOrientation','landscape'); hold on; % position [left bottom width height]
        position = 0.5;
        y_text_delta = 0.1;
        plot(my_output_fig7(5+2), position,'rd'); % Median
        plot([my_output_fig7(5+1); my_output_fig7(5+3)],(position)*[1;1],'r-') % NSE
        text(x_text, position+y_text_delta,'NSE interval (adj.)')
        position = 1;
        plot(my_output_fig7(1), position,'rd'); % ML - base case
        plot([my_output_fig7(2); my_output_fig7(3)],(position)*[1;1],'r-') % full interval with tau2_ML + trace(V)/k
        plot([my_output_fig7(4); my_output_fig7(5)],(position)*[1;1],'k|') % interval for mu
        text(x_text, position+y_text_delta,'Base case')
        ylim([0.2;1.3]);
        yticklabels({})
    end
    % to save to file go to the end of the script

    
    summaryStats_original = prctile(Y_i_original, [0 10 25 50 75 90 100],"Method","exact");
else

%% ML - case 4: equicorrelation - doesn't work, 
%%              replaced with autocorrelation

%fun4 = @(x)-phi_tau(Y_i, X, V_ec(V,-1/(k-1)+exp(x(2))), exp(x(1))); % equicorrelation: tau2 = exp(x(1)), rho = -1/(k-1)+exp(x(2))

fun4 = @(x)-phi_tau(Y_i, X, V_ac(V,-1+exp(x(2))), exp(x(1))); % autocorrelation: tau2 = exp(x(1)), rho = -1+exp(x(2))
if my_case<6
    x0 = [log(tau2_ML), log(0.9)];
else
    x0 = [log(0.01), log(0.5)];
end
[x_ML4_uc, fval] = fminunc(fun4,x0);
tau2_ML4_uc = exp(x_ML4_uc(1));
rho_ML4_uc = -1+exp(x_ML4_uc(2));
% use optimization with constraints
fun4 = @(x)-phi_tau(Y_i, X, V_ac(V,x(2)), x(1)); % tau2 = x(1), rho = x(2)
if my_case == 1
    x0 = [tau2_ML, 0.5];
elseif my_case < 6
    x0 = [tau2_ML, 0.9];
else
    x0 = [0.01, -0.1];
end
A = [];
b = [];
Aeq = [];
beq = [];
rho_lim = -1;% for equicorrelation -1/(k-1);
lb = [0, rho_lim];
ub = [Inf, 1];
[x_ML4, fval] = fmincon(fun4,x0,A,b,Aeq,beq,lb,ub);

tau2_ML4 = x_ML4(1);
rho_ML4 = x_ML4(2);
M_ML4 = inv(X' * inv(Omega(V_ac(V,rho_ML4),tau2_ML4)) * X) * (X' * inv(Omega(V_ec(V,rho_ML4),tau2_ML4)) * Y_i);
V_ML4_y = tau2_ML4 + trace(V)/k;
st_err_ML4 = sqrt(V_ML4_y);
V_ML4 = inv(X' * inv(Omega(V_ac(V,rho_ML4),tau2_ML4)) * X);
st_err_M_ML4 = sqrt(V_ML4);
if (my_case ~= 8) && (my_case ~= 6) 
    if (my_case == 1) || (my_case == 2) || (my_case == 4) || (my_case >= 9)
        position = -3;
    else
        position = -5;
    end
    plot(M_ML4, position,'rd'); % ML - case 4
    plot([M_ML4 - my_q*st_err_ML4; M_ML4 + my_q*st_err_ML4],(position)*[1;1],'r-') % full interval with tau2_ML + trace(V)/k
    %plot([M_ML4 - my_q*sqrt(trace(V)/k); M_ML4 + my_q*sqrt(trace(V)/k)],(-5)*[1;1],'kx') % interval with only trace(V)/k
    plot([M_ML4 - my_q*st_err_M_ML4; M_ML4 + my_q*st_err_M_ML4],(position)*[1;1],'k|') % interval for mu
    text(x_text, position+y_text_delta,'Correlated case')
end

% check likelihood as a function of rho for local optima
rho_set = [-1:0.01:1]';
phi_set =  rho_set * 0;
for j = 1:length(rho_set)
   phi_set(j) = phi_tau(Y_i, X, V_ac(V,rho_set(j)), tau2_ML4);
end
%figure;plot(rho_set, phi_set) % graph: concentrated likelihood as function of rho

%% extra case for revision: put autocorr structure to the other element

% use optimization with constraints
fun5 = @(x)-phi_tau1(Y_i, X, V, x(1), x(2)); % tau2 = x(1), rho = x(2)
x0 = [0.01, -0.1];
A = [];
b = [];
Aeq = [];
beq = [];
rho_lim = -1;% for equicorrelation -1/(k-1);
lb = [0, rho_lim];
ub = [Inf, 1];
[x_ML5, fval] = fmincon(fun5,x0,A,b,Aeq,beq,lb,ub);

tau2_ML5 = x_ML5(1);
rho_ML5 = x_ML5(2);
M_ML5 = inv(X' * inv(Omega1(V,tau2_ML5,rho_ML5)) * X) * (X' * inv(Omega1(V,tau2_ML5,rho_ML5)) * Y_i);
V_ML5_y = tau2_ML5 + trace(V)/k;
st_err_ML5 = sqrt(V_ML5_y);
V_ML5 = inv(X' * inv(Omega1(V,tau2_ML5,rho_ML5)) * X);
st_err_M_ML5 = sqrt(V_ML5);

position = -3.5;
plot(M_ML5, position,'rd'); % ML - extra case for revision 
plot([M_ML5 - my_q*st_err_ML5; M_ML5 + my_q*st_err_ML5],(position)*[1;1],'r-') % full interval with tau2_ML + trace(V)/k
plot([M_ML5 - my_q*st_err_M_ML5; M_ML5 + my_q*st_err_M_ML5],(position)*[1;1],'k|') % interval for mu
text(x_text+0.7, position+y_text_delta+0.1,'(a)')
text(x_text+0.7, position+y_text_delta-0.4,'(b)')

% check likelihood as a function of rho for local optima
rho_set = [-1:0.01:1]';
phi_set =  rho_set * 0;
for j = 1:length(rho_set)
   phi_set(j) = phi_tau1(Y_i, X, V, tau2_ML5, rho_set(j));
end
%figure; plot(rho_set, phi_set) % graph: concentrated likelihood as function of rho




%% collect all numbers together for Excel

my_output = [M_FE V_FE M_FE - my_q*sqrt(V_FE) M_FE + my_q*sqrt(V_FE); 
             M_RE V_RE M_RE - my_q*sqrt(V_RE) M_RE + my_q*sqrt(V_RE); 
             M_ML V_ML M_ML - my_q*st_err_M_ML M_ML + my_q*st_err_M_ML; 
             M_ML3 V_ML3 M_ML3 - my_q*st_err_M_ML3 M_ML3 + my_q*st_err_M_ML3;
             M_ML4 V_ML4 M_ML4 - my_q*st_err_M_ML4 M_ML4 + my_q*st_err_M_ML4
             M_ML5 V_ML5 M_ML5 - my_q*st_err_M_ML5 M_ML5 + my_q*st_err_M_ML5];
my_output_y = [V_ML_y M_ML - my_q*st_err_ML M_ML + my_q*st_err_ML; 
               V_ML3_y M_ML3 - my_q*st_err_ML3 M_ML3 + my_q*st_err_ML3;
               V_ML4_y M_ML4 - my_q*st_err_ML4 M_ML4 + my_q*st_err_ML4
               V_ML5_y M_ML5 - my_q*st_err_ML5 M_ML5 + my_q*st_err_ML5];
my_output_table1 = [k sqrt(tau2_ML) sqrt(trace(V)/k) st_err_ML (trace(V)/k)/V_ML_y];

my_output_table2 = [sqrt(tau2_ML3) sqrt(trace(sigma2_eps*V)/k) st_err_ML3 (trace(sigma2_eps*V)/k)/V_ML3_y];

my_output_table3a = [sqrt(tau2_ML4) sqrt(trace(V)/k) st_err_ML4 (trace(V)/k)/V_ML4_y rho_ML4];

my_output_table3b = [sqrt(tau2_ML5) sqrt(trace(V)/k) st_err_ML5 (trace(V)/k)/V_ML5_y rho_ML5];
end

%% print into file
print_flag = 0;
if (print_flag == 1) 
    path_name = ['/Users/av/Dropbox (Sydney Uni)/Research/2023/IPCC/Draft_meta/figures/']; % mac path name
    %path_name = ['C:\Users\avasnev\Dropbox (Sydney Uni)\Research\2023\IPCC\Draft_meta\figures\']; % PC path name
    %path_name = []; % online
    if my_case < 12
        output_name = [path_name 'Figure_meta_case_' num2str(my_case) '_revision.eps'];   
    else
        output_name = [path_name 'Figure_meta_case_' num2str(my_case) '_RT-H' num2str(my_hypothesis) '_Stage' num2str(my_stage) '.eps']; % for figure with cases
        %output_name = [path_name 'Figure_meta_case_' num2str(my_case) '_RT-H' num2str(my_hypothesis) '_Stage' num2str(my_stage) '_2intervals.eps']; % for figure with only two intervals
    end
    print(output_name,'-dpsc2');
    %print(output_name,'-depsc'); % formats '-dpsc2' - I used before; '-dpdf' - doesn't have bounding box
    %exportgraphics(gcf,output_name,'BackgroundColor','none','ContentType','vector')
end 

%% functions

function err = err(y, X, V, tau2)
   beta_hat = inv(X' * inv(Omega(V,tau2)) * X) * (X' * inv(Omega(V,tau2)) * y);
   err = y - X * beta_hat;
end

function err1 = err1(y, X, V, tau2, rho) % extra case for revision
   beta_hat = inv(X' * inv(Omega1(V,tau2,rho)) * X) * (X' * inv(Omega1(V,tau2,rho)) * y);
   err1 = y - X * beta_hat;
end

function Omega = Omega(V,tau2)
   Omega = tau2 * eye(size(V,1)) + V;
end

function Omega1 = Omega1(V,tau2,rho) % extra case for revision
   %Omega = tau2 * eye(size(V,1)) + V; % structure used before the revision
   k = size(V,1);
   P_ac = eye(k); % autocorrelation structure
   for j = 1:k-1
       P_ac = P_ac + rho^j * (diag(ones(k-j,1),j) + diag(ones(k-j,1),-j));
   end
   P = P_ac; % test autocorrelation structure
   Omega1 = tau2 * P + V;
end

function phi_tau = phi_tau(y, X, V, tau2)
   tmp = err(y, X, V, tau2)' * inv(Omega(V, tau2)) * err(y, X, V, tau2);
   phi_tau = -(log(det(Omega(V, tau2))) + tmp);
   % phi_tau = -(tmp); % try least squares
end

function phi_tau1 = phi_tau1(y, X, V, tau2, rho) % extra case for revision
   tmp = err1(y, X, V, tau2, rho)' * inv(Omega1(V, tau2, rho)) * err1(y, X, V, tau2, rho);
   phi_tau1 = -(log(det(Omega1(V, tau2, rho))) + tmp);
end

function V_ec = V_ec(V,rho)
   %rho = 0.2; % manual override to check the fixed value of rho 
   k = size(V,1);
   P = (ones(k,k) - eye(k)) * rho + eye(k); % equicorrelation
   tmp = sqrtm(V);
   V_ec = tmp * P * tmp;
end

function V_ac = V_ac(V,rho)
   %rho = 0.2; % manual override to check the fixed value of rho 
   k = size(V,1);
   %P = (ones(k,k) - eye(k)) * rho + eye(k); % equicorrelation
   P_ac = eye(k); % autocorrelation structure
   for j = 1:k-1
       P_ac = P_ac + rho^j * (diag(ones(k-j,1),j) + diag(ones(k-j,1),-j));
   end
   P = P_ac; % test autocorrelation structure
   tmp = sqrtm(V);
   V_ac = tmp * P * tmp;
end

function remove_idx = remove_LCS(dataTable)
remove_idx = [];
for my_hypothesis = 1:6
    for my_stage = 1:4
        % select slice of the data for one hypothesis and one stage
        my_selection = (dataTable.rt_hypothesis == my_hypothesis)&(dataTable.stage == my_stage);
        dataTable_filtered = dataTable(my_selection,:);
        Y_i = dataTable_filtered.estimate;

        remove_pct = 0.05; % percentage to remove
        n_extr = round(remove_pct * length(Y_i)); % number of extreme observations to remove
        % Sort and get original indices
        [~, sorted_idx] = sort(Y_i);
        % Indices of 5% smallest
        idx_smallest_5pct = sorted_idx(1:n_extr);
        % Indices of 5% largest
        idx_largest_5pct = sorted_idx(end-n_extr+1:end);
        
        % also remove outliers in st deviation
        sigma_i = dataTable_filtered.standard_error;
        [~, sorted_idx_sigma] = sort(sigma_i);
        idx_smallest_5pct_sigma = sorted_idx_sigma(1:n_extr);
        idx_largest_5pct_sigma = sorted_idx_sigma(end-n_extr+1:end);
        
        remove_idx = union(remove_idx, idx_smallest_5pct); % combine both sets
        remove_idx = union(remove_idx, idx_largest_5pct);  
        remove_idx = union(remove_idx, idx_largest_5pct_sigma);  % add outliers in st dev
        remove_idx = union(remove_idx, idx_smallest_5pct_sigma);
    end
end
end