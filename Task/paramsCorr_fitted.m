%% Correlation between fitted parameters on subjects
detRGB = [0, 0.4470, 0.7410];
infRGB = [0.8500, 0.3250, 0.0980];
selRGB = [0.4660, 0.6740, 0.1880];

n_subjects  = 30;
excluded    = [14 20 21 22 27]; % excluded subject numbers

load('glazeOrigParamsLog.mat');
load('infNoiseParamsLog');
load('selNoiseParamsLog');

h_det = [];
h_inf = [];
h_sel = [];

infN_inf = [];
infN_sel = [];

selN = [];

for i = 1:n_subjects
    % skip excluded subject numbers
    if ismember(i,excluded)
        continue;
    end
    
    h_det = vertcat(h_det, glazeOrigParamsLog(i).past);
    h_inf = vertcat(h_inf, infNoiseParamsLog(i).past(1));
    h_sel = vertcat(h_sel, selNoiseParamsLog(i).past(1));
    
    infN_inf = vertcat(infN_inf, infNoiseParamsLog(i).past(2));
    infN_sel = vertcat(infN_sel, selNoiseParamsLog(i).past(2));
    
    selN     = vertcat(selN, selNoiseParamsLog(i).past(3));
end

% Compare hazard rate parameters
[rsqrd_det_inf_H(1) rsqrd_det_inf_H(2)]  = corr(h_det, h_inf); %pearson correlation
[rsqrd_det_sel_H(1) rsqrd_det_sel_H(2)]  = corr(h_det, h_sel); %pearson correlation
[rsqrd_inf_sel_H(1) rsqrd_inf_sel_H(2)]  = corr(h_inf, h_sel); %pearson correlation
% Get coeffients for line of best fit
det_inf_H_coeffs = polyfit(h_det, h_inf,1); %linear coeffs
det_sel_H_coeffs = polyfit(h_det, h_sel,1); %linear coeffs
inf_sel_H_coeffs = polyfit(h_inf, h_sel,1); %linear coeffs

% Compare inference noise parameters
[rsqrt_inf_sel_infN, infN_pval] = corr(infN_inf, infN_sel); %pearson correlation
infN_coeffs                     = polyfit(infN_inf, infN_sel,1); %linear coeffs

% Correlation of the selection noise parameter with change in inference noise
% parameter
diff_infN           = infN_sel - infN_inf;    %difference in the infNoise parameter
rsqrd_sel_inf_selN  = corr(diff_infN, selN) %pearson correlation
selN_coeffs         = polyfit(diff_infN, selN,1); %linear coeffs


%% Plotting
close all;
figure(1);
%title('Correlation of fitted H parameters of models on subjects in postdiction condition');
hold on;
det_inf_H = scatter(h_det, h_inf);
det_sel_H = scatter(h_det, h_sel,'x');
inf_sel_H = scatter(h_inf, h_sel,'d');
colorOrder = get(gca,'ColorOrder');
plot([0:1], polyval(det_inf_H_coeffs,[0:1]),'Color',colorOrder(1,:));
plot([0:1], polyval(det_sel_H_coeffs,[0:1]),'Color',colorOrder(2,:));
plot([0:1], polyval(inf_sel_H_coeffs,[0:1]),'Color',colorOrder(3,:));
plot([0:1], [0:1],'LineStyle',':','Color',[0 0 0])
text(.43,.14, ['r^2 (Pearson) = ' num2str(rsqrd_det_inf_H)],'Color',colorOrder(1,:));
text(.43,.12, ['r^2 (Pearson) = ' num2str(rsqrd_det_sel_H)],'Color',colorOrder(2,:));
text(.43,.10, ['r^2 (Pearson) = ' num2str(rsqrd_inf_sel_H)],'Color',colorOrder(3,:));
legend([det_inf_H det_sel_H inf_sel_H], {'Inf. Noise to Determinstic', 'Sel. Noise to Deterministic', 'Sel. Noise to Inf. Noise'})
xlabel('H (relative baseline)');
ylabel('H (comparison)');
xlim([0 .6]);
ylim([0 .6]);
hold off;

figure(2);
%title({'Correlation of fitted inference noise parameters (\sigma) of models on subjects in postdiction condition', ' '});
hold on;
inf_sel_infN = scatter(infN_inf, infN_sel, 'MarkerEdgeColor', infRGB);
plot([0:2], [0:2],'LineStyle',':','Color',[0 0 0]);
plot([0:2], polyval(infN_coeffs,[0:2]),'Color',infRGB);
text(.9,.6, ['r^2 (Pearson) = ' num2str(rsqrt_inf_sel_infN)],'Color',infRGB);
%legend([inf_sel_infN], {'Sel. Noise to Inf. Noise'});
xlabel('\sigma (model with inference noise only)');
ylabel('\sigma (model with inference noise and selection noise)');
xlim([0 1.2]);
ylim([0 1.2]);
hold off;

figure(3);
%title({'Correlation of fitted selection noise parameters (\eta) over the change in inference noise parameters /non subjects in postdiction condition', ' '});
hold on;
scatter(diff_infN, selN,'MarkerEdgeColor',selRGB);
plot([-1:1], polyval(selN_coeffs,[-1:1]),'Color',selRGB);
text(-.5,1, ['r^2 (Pearson) = ' num2str(rsqrd_sel_inf_selN)],'Color',selRGB);
xlabel('Change in inference noise parameter (\sigma(mixed-noise) - \sigma only)');
ylabel('\eta (Selection Noise)');
xlim([-.7 .1]);
hold off;

