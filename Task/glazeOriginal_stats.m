% Obtain statistics from the parameters obtained with BADS
%
%   Name:   glazeOriginal_stats
%   Type:   Script
% Author:   Jun Seok Lee
%   Date:   March 2019
%   Note:   Adapted for TURFU experiment
%
% This script calls upon 2 files:
%   1. glazeOriginal.m
%   2. glazeOriginal_paramSearchBADS.m

global globCond;
global globSubj;
global blockFilter;

% Specify whether we are trying to recover a model
modelRec = false;

% Set the blocking condition on which to estimate parameters 
blockFilter = 'direction'; % or set to 'volatility'

cond = [1 2];
nSubj = 29;
maxIters = 5;
excluded = [6 8 14 18 20 21];

backCondHazRts = zeros(1,maxIters);
forwCondHazRts = zeros(1,maxIters);
subjHazRts = cell(1,24);

if strcmpi(blockFilter, 'direction')
    cond1txt = 'postdictive inference';
    cond2txt = 'predictive inference';
elseif strcmpi(blockFilter, 'volatility')
    cond1txt = 'low volatility';
    cond2txt = 'high volatility';
else
    err("Choose appropriate blocking variable: 'condition' or 'volatility'");
end
    

for isubj = 1:nSubj
    
    if ismember(isubj, excluded) 
        continue;   % Don't analyze subjects that are in the exclusion group
    end  
    
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
            run glazeOriginal_paramSearchBADS;
            if i == 1
                backCondHazRts(j) = optiParams;
            else
                forwCondHazRts(j) = optiParams;
            end
        end
    end
    subjHazRts{isubj} = {backCondHazRts forwCondHazRts};
end

save('subjPriorHazRts_glazeOriginal', 'subjHazRts');

%% Model simulation
% This takes the parameters obtained from the above script and simulates the behavior
% of a subject corresponding to those parameters

%% Model recovery
% This should essentially run the initial script above to see if we get back the same
% parameters back. 

%% Intra-subject "significance testing" for differences in mean (single subject analysis)

tailedStats = zeros(1,24); 
Stats = zeros(1,24);       

for isubj = 1:nSubj
% H1: Do subjects perceive lower volatility in the agent condition?
    [p, h, stats] = ranksum(subjHazRts{isubj}{1}, subjHazRts{isubj}{2}, 'Tail','right'); %wilcoxon rank sum
    tailedStats(isubj) = h;
end
sigrate = numel(tailedStats(tailedStats ==1))/24;
disp(['Stat. sig. for ' num2str(sigrate) ' of subjects being MORE reliant on prior in agent condition.']);

%% Inter-subject difference significance testing group level (group level analysis)
diff = [];
for isubj = 1:nSubj
    diff = horzcat(diff, mean(subjHazRts{isubj}{2}) - mean(subjHazRts{isubj}{1})); %wilcoxon rank sum
end
[p, h, stats] = ranksum(diff,zeros(1,24),'Tail','left')
