% Model simulator script for the Forward Glaze (2015) model + Inference Noise
% 
%   Name:   infNoiseGlazeForward_sim.m
% Author:   Jun Seok Lee
%   Date:   April 2018
%   Type:   Script
%   Note:   For the TURFU experiment


% need to load in fitted parameters (H, sigma), experiment values (seqllrs, nsamples)

% Load parameters obtained from model fitting
load('paramsLog','paramsLog'); % this could be changed based on which parameter set is being used

% Static variables (set and forget)
n_subjects  = 30;   % set to number of subjects we are analyzing
excluded    = [21]; % excluded subject numbers

exclude     = true; % set to true to exclude shitty results
if exclude
    excluded = [14 20 21 22 27];
end
subjects = [1:n_subjects];

% Set up data structure for the generated responses
simBehaviorForward = struct;

modelLL = zeros(30,4); % structure to hold model likelihoods to subject data

for isubj = subjects
    % skip excluded subjects
    if ismember(isubj, excluded)
        continue;
    end
    
    % assign shorter variable names for convenience
    haz     = paramsLog(isubj).past(1);
    sigm    = paramsLog(isubj).past(2);
    
    simBehaviorForward(isubj).pastParams    = paramsLog(isubj).past;
    simBehaviorForward(isubj).rslt          = struct;
    
    % load experimental behavior from subjects
    filename = dir(sprintf('Data/TURFU_S%02d_*.mat', isubj));
    filename = filename.name;
    load(sprintf('Data/%s',filename));
    
    % Find blocks corresponding to the future condition & not training
    blocks = find([expe.blck.taskid] == 2 & [expe.blck.condtn] ~= 3);
    
    iblock = 1; % iterative ctr on blocks
    for block = blocks
        seqllrs     = expe.blck(block).seqllr;
        nsamples    = [];
        % loop through the trials in current block to get number of samples
        for itrial = 1:72 
            nsamples = horzcat(nsamples, numel(expe.blck(block).seqtlt{itrial}));
        end
        
        % simulate behavior in future blocks based on parameters and experimental conditions
        simBehaviorForward(isubj).rslt(block).resp = infNoiseGlazeForward_simulator(haz, sigm, seqllrs, nsamples);
        % unnecessary code below, but for aesthetics in the simBehaviorForward structure
        if ~ismember(9,blocks)
            simBehaviorForward(isubj).rslt(9).resp  = [];
            simBehaviorForward(isubj).rslt(10).resp = [];
        end
        
        corTrls = simBehaviorForward(isubj).rslt(block).resp(2:73)-expe.rslt(block).resp(2:73); % correct generations
        modelLL(isubj,iblock) = log(numel(corTrls(corTrls==0))/72); % take log of likelihood of model to behavior on block
        
        % DEBUG: code for visualizing the generated responses to subjects' %%
        if false % set to true for DEBUG
            figure;
            for jtrial = 2:73
                xlim([0 73]);
                ylim([.5 2.5]);
                scatter([1:72], simBehaviorForward(isubj).rslt(block).resp(2:73));
                hold on;
                scatter([1:72], expe.rslt(block).resp(2:73)-.1);
                hold off;
            end
            waitforbuttonpress;
            close all;
        end
        %%%%%%%%%%%
        iblock = iblock+1;
    end
    
end

% Calculate total likelihood that the data came from the model
totalLL = 0;
for isubj = 1:numel(modelLL(:,1))
    if modelLL(isubj,1) == 0
        continue;
    end
    
    totalLL = totalLL + sum(modelLL(isubj,:)); 
end
totalLL = totalLL.*-1;

