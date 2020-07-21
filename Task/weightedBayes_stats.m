% Obtain statistics from the parameters obtained with BADS
%
% Type:     Script
% Author:   Jun Seok Lee
% Date:     February 2019
%
% This script calls upon 2 files:
%   1. weightedBayes.m
%   2. weightedBayes_paramSearchBADS.m

global globCond;
global globSubj;
global blockFilter;

% Specify whether we are trying to recover a model
modelRecovery = false;

% Set the blocking condition on which to estimate parameters 
blockFilter = 'condition'; % or set to 'volatility'

cond = [1 2];
nSubj = 24;
maxIters = 10;

backCondWts = zeros(1,maxIters);
forwCondWts = zeros(1,maxIters);
subjWts = cell(1,24);

if blockFilter == 'condition'
    cond1txt = 'observer/backward inference'
    cond2txt = 'agent/forward inference'
elseif blockFilter == 'volatility'
    cond1txt = 'low volatility';
    cond2txt = 'high volatility';
else
    err("Choose appropriate blocking variable: 'condition' or 'volatility'");
end
    

for isubj = 1:nSubj
    globSubj = isubj;
    for i = cond
        globCond = i;
        % Print block condition for verification before program start
        if cond == 1 
            disp(['Analyzing ' cond1txt ' blocks...']);
        else
            disp(['Analyzing ' cond2txt ' blocks...']);
        end
        for j = 1:maxIters
            run weightedBayes_paramSearchBADS;
            if i == 1
                backCondWts(j) = optiParams;
            else
                forwCondWts(j) = optiParams;
            end
        end
    end
    subjWts{isubj} = {backCondWts forwCondWts};
end

save('subjPriorWts', 'subjWts');

%% Intra-subject "significance testing" for differences in mean (single subject analysis)

tailedStats = zeros(1,24); 
Stats = zeros(1,24);       

for isubj = 1:nSubj
% H1: Are subjects more reliant on the prior in the agent condition?
    [p, h, stats] = ranksum(subjWts{isubj}{1}, subjWts{isubj}{2}, 'Tail','left'); %wilcoxon rank sum
    tailedStats(isubj) = h;
end
sigrate = numel(tailedStats(tailedStats ==1))/24;
disp(['Stat. sig. for ' num2str(sigrate) ' of subjects being MORE reliant on prior in agent condition.']);

%% Inter-subject difference significance testing group level (group level analysis)
diff = [];
for isubj = 1:nSubj
    diff = horzcat(diff, mean(subjWts{isubj}{2}) - mean(subjWts{isubj}{1})); %wilcoxon rank sum
end
[p, h, stats] = ranksum(diff,zeros(1,24),'Tail','right')
