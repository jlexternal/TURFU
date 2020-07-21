% Obtain Data|Model log-likelihood 
%
%   Name:   modelLL_get.m
%   Type:   Script
% Author:   Jun Seok Lee
%   Date:   May 2019
% 
% Description:  Calculates the actual log-likelihood of data|model for the stochastic
%               models.
%               Use this code ONLY after the parameters have been found using the
%               LogSumExp workaround with the BADS algo.

global blockFilter;
global globCond;
global globParticleCount;

blockFilter = 'direction';
globParticleCount = 5000;
globCond = 2;       % 1 past, 2 future
if globCond == 2
    futurecond = true;  % true for future condition
end

n_subjects = 30;
excluded = [14 20 21 22 27];

load('infNoiseParamsLog');
load('selNoiseParamsLog');
modelLLs = [];

isubjCtr = 1;
for isubj = 1:n_subjects
    if ismember(isubj, excluded)
       continue; 
    end
    filename = dir(sprintf('Data/TURFU_S%02d_*.mat', isubj));
    filename = filename.name;
    load(sprintf('Data/%s',filename));
    
    params = infNoiseParamsLog(isubj).past;
    modelLLs(isubjCtr,1) = infNoiseGlaze_PF_recov(params,true,futurecond);
    
    params = selNoiseParamsLog(isubj).past;
    modelLLs(isubjCtr,2) = selNoiseGlaze_PF_recov(params,true,futurecond); 
    
    isubjCtr = isubjCtr + 1;
end

% save manually here 
% using save('modelLLs','modelLLs');

%% Calculate BIC

%load('modelLLs');
bic = @(n, k, ll) (log(n)*k) - (2*ll);
k       = [2 3];    % number of parameters for inf and sel noise models
n       = 288;      % number of trials per inferential direction condition (postdictive)
bic_inf = 0;
bic_sel = 0;
for i = 1:size(modelLLs,1)
    bic_inf = bic_inf + bic(n, k(1), -modelLLs(i,1));
    bic_sel = bic_sel + bic(n, k(2), -modelLLs(i,2));
end

disp(['BIC (Inf. Noise) = ' num2str(bic_inf)]);
disp(['BIC (Inf. Noise + Sel. Noise) = ' num2str(bic_sel)]);
