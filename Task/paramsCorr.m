%% Fitted vs Recovered Parameters (on actual subject data)
% This takes the fitted parameters on subjects and compares them to the recovered
% parameters. Can be run on either the deterministic or the noisy model.

n_subjects  = 30;
excluded    = [14 20 21 22 27]; % excluded subject numbers

hRate_past = [];
hRate_past_recd = [];
sigm_past = [];
sigm_past_recd = [];

model = 'glaze'; % 'glaze' or 'infNoise' or 'selNoise'

for i = 1:n_subjects
    % skip excluded subject numbers
    if ismember(i,excluded)
        continue;
    end
    
    if strcmpi(model, 'infNoise')
    hRate_past      = vertcat(hRate_past, paramsLog(i).past(1));
    hRate_past_recd = vertcat(hRate_past_recd, paramsLog(i).pastRecd(1));
    sigm_past       = vertcat(sigm_past, paramsLog(i).past(2));
    sigm_past_recd  = vertcat(sigm_past_recd, paramsLog(i).pastRecd(2));
    elseif strcmpi (model, 'glaze')
    hRate_past      = vertcat(hRate_past, glazeOrigParamsLog(i).past(1));
    hRate_past_recd = vertcat(hRate_past_recd, glazeOrigParamsLog(i).pastRecd(1));
    sigm_past       = vertcat(sigm_past, glazeOrigParamsLog(i).past(1));
    sigm_past_recd  = vertcat(sigm_past_recd, glazeOrigParamsLog(i).pastRecd(1));
    end
   
    
end

if strcmpi(model, 'infNoise')
    titletext = 'Hazard rate (H) parameters (Recovered vs. Fitted) on Glaze (2015) + Inference Noise Model';
    hRate_past_sto = hRate_past;
elseif strcmpi (model, 'glaze')
    titletext = 'Hazard rate (H) parameters (Recovered vs. Fitted) on Glaze (2015) Deterministic Model';
    hRate_past_det = hRate_past;
end

clf;
figure(1);
hold on;
title(titletext);
xlabel('H from fitting models to subject behaviors'' past condition');
ylabel('Recovered H');
xlim([0 .6]);
ylim([0 .6]);
scatter(hRate_past, hRate_past_recd);
coeffs = polyfit(hRate_past, hRate_past_recd,1);
rsqrd = fitlm(hRate_past, hRate_past_recd);
fittedX = linspace(0, .6, 200);
fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY,'Color', 'r','LineWidth',.6);	% regression line
plot(fittedX, fittedX,':','Color','black');         % pure diagonal
text(fittedX(150)+.06, fittedY(150)+.01, ['r^2 = ' num2str(rsqrd.Rsquared.Ordinary)],'Color','r');
hold off;


if strcmpi(model, 'infNoise')
    figure(2);
    hold on;
    title('Inference Noise (\sigma) parameter (Recovered vs. Fitted to subjects'' past condition)');
    xlabel('\sigma from fitting models to subjects'' past condition');
    ylabel('Recovered \sigma');
    xlim([0 1]);
    ylim([0 1]);
    scatter(sigm_past, sigm_past_recd);
    coeffs = polyfit(sigm_past, sigm_past_recd,1);
    rsqrd = fitlm(sigm_past, sigm_past_recd);
    fittedX = linspace(0, 1, 200);
    fittedY = polyval(coeffs, fittedX);
    plot(fittedX, fittedY,'Color', 'r','LineWidth',.6); % regression line
    plot(fittedX, fittedX,':','Color','black');         % pure diagonal
    text(fittedX(150)+.06, fittedY(150)+.01, ['r^2 = ' num2str(rsqrd.Rsquared.Ordinary)],'Color','r');
    hold off;
end

%% Cross-model (Deterministic -> Noisy) Comparison on Hazard Rate
% This show the relationship between the hazard rate parameter on the model with and
% without noise.

figure(3);
hold on;
title('Hazard Rates (H) from the Stochastic model vs. Deterministic (past condition)');
xlabel('H fitted to subjects with Deterministic (purely Glaze) Model');
ylabel('H fitted to subjects with Stochastic (Glaze + Inf. Noise) Model');
xlim([0 1]);
ylim([0 1]);
scatter(hRate_past_det, hRate_past_sto);
coeffs = polyfit(hRate_past_det, hRate_past_sto,1);
rsqrd  = fitlm(hRate_past_det, hRate_past_sto);  %pearson correlation
spcorr = corr(hRate_past_det, hRate_past_sto,'Type','Spearman');
fittedX = linspace(0, 1, 200);
fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY,'Color', 'r','LineWidth',.6); % regression line
plot(fittedX, fittedX,':','Color','black');         % pure diagonal
text(fittedX(150)-.05, fittedY(150)-.09, ['r^2 (Pearson) = ' num2str(rsqrd.Rsquared.Ordinary)],'Color','r');
text(fittedX(150)-.05, fittedY(150)-.12, ['r^2 (Spearman) = ' num2str(spcorr)],'Color','b');
hold off;

%% Artificial values recovery
% This is for analysis on the recovered parameters from determined values run on several
% experimental conditions.

n_blocks = numel(simBehavior(1).expe);
allH = [];
hIndex = 1;
for hval = [simBehavior.H_val]
    %calculations
    allH = horzcat(allH, [ones(1,n_blocks).*hval;simBehavior(hIndex).H_val_recd]);
    
    avgH = mean(simBehavior(hIndex).H_val_recd);
    stdH = std(simBehavior(hIndex).H_val_recd);
    %plotting
    scatter(ones(1,n_blocks).*hval, simBehavior(hIndex).H_val_recd,'filled','MarkerFaceColor',[255./255 206./255 206./255],...
            'LineWidth', 2);
    errorbar(hval, avgH, stdH, 'Color',[255./255 88./255 88./255],'LineWidth',1,'CapSize',0);
    scatter(hval, avgH, 'r');
    hold on;
    hIndex = hIndex + 1;
end
% calculations
coeffs_H_artif = polyfit(allH(1,:), allH(2,:),1);
rsqrd_H_artif  = fitlm(allH(1,:), allH(2,:));  %pearson correlation
%plotting
title({'Recovered vs Generative parameters (H) on the Deterministic Glaze (2015)', 'model run on 8 experimental conditions'});
xlabel('Generative H');
ylabel('Recovered H');
plot([0 .7], [0 .7],':', 'Color',[0 0 0],'LineWidth',1.5);
fittedX = linspace(0, .7, 20);
fittedY = polyval(coeffs_H_artif, fittedX);
plot(fittedX, fittedY,'Color', 'r','LineWidth',.6); % regression line
text(fittedX(10)-.05, fittedY(10)-.09, ['r^2 (Pearson) = ' num2str(rsqrd.Rsquared.Ordinary)],'Color','r');
hold off;



