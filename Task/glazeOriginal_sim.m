% Model simulator script for the Glaze (2015) model 
%
%   Name:   glazeOriginal_sim.m
% Author:   Jun Seok Lee
%   Date:   April 2018
%   Type:   Script
%   Note:   Adapted for TURFU experiment

clearvars -except infNoiseMeanRespRate infNoiseInclRespRate inclRespRate meanRespRate;

% Global variables
global blockFilter;

% Load parameters obtained from model fitting
load('glazeOrigParamsLog','glazeOrigParamsLog');

% BADS settings
addpath(genpath('bads-master')); % load the BADS files into working directory

% Search Parameters
%   Given that sigma and hazard rate are our parameters, the vectors below hold
%   values in the order: [hazRate]
startPt = [rand()]; % vector corresponding to values of parameters to be estimated
                    %   Side note: should probably set this to various functions on the
                    %    parameter map to see if the optimum is found around the same
                    %    values on the noisy surface
lBound  = [0];      % HARD lower bound of parameter values
uBound  = [1];      % HARD upper bound of parameter values
pLBound = [.01];    % Plausible lower bound of parameter values (optional)
pUBound = [.99];    % Plausible upper bound of parameter values (optional)

% Changing parameter tolerance options on BADS
options = [];
options.TolMesh = 0.01;

% Static variables (set and forget)
n_subjects  = 30;   % set to number of subjects we are analyzing
excluded    = [21]; % excluded subject numbers

if n_subjects ~= numel(glazeOrigParamsLog) % error catching
    error('Number of subjects is not equal to the number of entries in the parameters log.');
end

% Dynamic variables (CHANGE VALUES HERE BASED ON ANALYSIS)
blockFilter     = 'direction';  % 'direction' or 'volatility'
paramRecovery   = false;        % set to true if trying to recover parameters
exclude         = true;         % set to true to exclude shitty results
specSubjects    = false;        % set to true if analyzing only certain subjects=
condition      = [1];       % set to 1 if for 1st condition, 2 for 2nd, 1:2 for both
disp(['Running the ' num2str(condition) 'st/nd ' blockFilter ' condition blocks...']);

if exclude
    excluded = [14 20 21 22 27];
end
if specSubjects
    subjects = [10];          % specify subjects here
else
    subjects = [1:n_subjects];
end

% Data structures to hold simulated behavior
% In this deterministic model, each subject will have 2 fields nested under the rslt
% field: 1) resp,   for the actual responses
%        2) belief, for the posterior belief LLR leading to the response
simBehavior = struct;

% stuff for model fit to data, specifically p(choice2)|seqllr
glazeorigchoice1     = zeros(30,8);
glazeorigchoicetot   = zeros(30,8);

modelLL = zeros(30,4); % structure to hold model likelihoods to subject data

%% Simulate behavior using parameters from subjects

for isubj = subjects
    % skip excluded subjects
    if ismember(isubj, excluded)
        continue;
    end
    % load experimental behavior from subjects
    filename = dir(sprintf('Data/TURFU_S%02d_*.mat', isubj));
    filename = filename.name;
    load(sprintf('Data/%s',filename));
    
    % identify condition that parameters were found in
    for cond = condition
        
        % Pseudocode:
        % 1. Find the 4 blocks and their seqllrs that the subject saw 
        % 2. Simulate behavior based on those blocks and the fitted parameters
        % 3. Run the particle filter on the behavior and the experimental settings
        % 4. Compare the parameters
        
        
        % Reference for condition and fitted parameter
        if strcmpi(blockFilter, 'direction')
            blocks = find([expe.blck.taskid] == cond);	% taskid refers to inference direction in exp.
            haz     = glazeOrigParamsLog(isubj).past;
        elseif strcmpi(blockFilter, 'volatility')
            blocks = find([expe.blck.condtn] == cond);	% condtn refers to volatility of block
            if cond == 1
                haz     = glazeOrigParamsLog(isubj).low;
            else
                haz     = glazeOrigParamsLog(isubj).high;
            end
        elseif strcmpi(blockFilter, 'all')
            blocks = [3:10];    % when analyzing all non-training blocks
        else
            error("Choose appropriate blocking variable: 'condition' (inference direction) or 'volatility'")
        end
        blocks(find(blocks==1)) = [];   % disregard training blocks
        blocks(find(blocks==2)) = [];   % disregard training blocks
        
        iblock = 1;
        % loop through the blocks corresponding to current condition (absolute index; do not increment)
        for block = blocks
            seqllrs     = expe.blck(block).seqllr;
            nsamples    = [];
            
            % loop through the trials in current block to get number of samples
            for itrial = 1:72 
                nsamples = horzcat(nsamples, numel(expe.blck(block).seqtlt{itrial}));
            end
            
            % simulate behavior based on params for block and store into local memory
            if strcmpi(blockFilter, 'direction') && condition == 1
                [simBehavior(isubj).rslt(block).resp, simBehavior(isubj).rslt(block).belief] = glazeOriginal_simulator(haz, seqllrs, 0);
                                                                                               % parameters(haz, seqllrs, futurecond)
            elseif strcmpi(blockFilter, 'direction') && condition == 2
                [simBehavior(isubj).rslt(block).resp, simBehavior(isubj).rslt(block).belief] = glazeOriginal_simulator(haz, seqllrs, 1);
            end
            
            if paramRecovery
                % RISKY CODE (need to verify that this does not change the save file)
                % replace human responses with simulated responses directly into expe variable
                expe.rslt(block).resp = simBehavior(isubj).rslt(block).resp; 
            end
            
            % Model fit to data stuff
            for itrial = 1:72
                tempResp = simBehavior(isubj).rslt(block).resp(itrial+1);
                tempChoicetot = glazeorigchoicetot(isubj,binIndex(expe.blck(block).seqllr(itrial)));
                glazeorigchoicetot(isubj,binIndex(expe.blck(block).seqllr(itrial))) = tempChoicetot + 1;
                if tempResp == 1
                    tempChoice1 = glazeorigchoice1(isubj,binIndex(expe.blck(block).seqllr(itrial)));
                    glazeorigchoice1(isubj,binIndex(expe.blck(block).seqllr(itrial))) = tempChoice1 + 1;
                end
            end
            
            % model log-likelihood calculation
            corTrls = simBehavior(isubj).rslt(block).resp(2:73)-expe.rslt(block).resp(2:73); % correct generations
            modelL(isubj,iblock) = numel(corTrls(corTrls==0))/72; % take likelihood of model to behavior on block for isubj
        
            iblock = iblock + 1;
        end
        
        if paramRecovery
            % Initializing values for BADS visuals
            txt = 1;
            x1 = [];
            clf;
            % run BADS
            [optiParams, fVal] = bads(@glazeOriginal, startPt, lBound, uBound, pLBound, pUBound, [], options)

            % Register the recovered parameters into the original paramsLog structure
            if strcmpi(blockFilter, 'direction')
                if cond == 1
                    glazeOrigParamsLog(isubj).pastRecd   = optiParams;
                else
                    glazeOrigParamsLog(isubj).futureRecd = optiParams;
                end
            elseif strcmpi(blockFilter, 'volatility')
                if cond == 1
                    glazeOrigParamsLog(isubj).lowRecd    = optiParams;
                else
                    glazeOrigParamsLog(isubj).highRecd   = optiParams;
                end
            end
        end
    end
end

% save simulated choices and beliefs to file
dateStr = datestr(datenum(date),'ddmmyyyy');
if strcmpi(blockFilter,'direction')
    if condition == 1
        save(['glazeOriginal_sim_post_struct_' dateStr], 'simBehavior');
    elseif condition == 2
        save(['glazeOriginal_sim_pred_struct_' dateStr], 'simBehavior');
    end
end

% Calculate total likelihood that the data came from the model
totalL = [];
for isubj = 1:numel(modelL(:,1))
    if modelL(isubj,1) == 0
        continue;
    end
    % list the mean likelihood across the blocks for each subject
    totalL = vertcat(totalL, sum(modelL(isubj,:))/numel(modelL(isubj,:))); 
end
totalL  = mean(totalL); % take mean likelihood across the subjects
totalLL = log(totalL);  % take the log of the mean likelihood

%% Convert to percentages and plot

glazeOrigRespRate = glazeorigchoice1./glazeorigchoicetot;
exclRows = isnan(glazeOrigRespRate(:,1));
glazeOrigInclRespRate = glazeOrigRespRate(~exclRows,:);
glazeOrigMeanRespRate = mean(glazeOrigInclRespRate);

scatter([1:8], glazeOrigMeanRespRate);
hold on;
errorbar([1:8], glazeOrigMeanRespRate, std(glazeOrigInclRespRate));
hold off;

figure(2);
hold on;
for itest = 1:25
    scatter([1:8],glazeOrigInclRespRate(itest,:));
end
hold off;

%% Local Functions
function binNum = binIndex(signllr)
    % This function assigns the bin number depending on the llr value signed by
    %   the previous subject choice

    % data structure: binStruct: [(-inf,-3), [-3:-2), [-2:-1), [-1:0), [0:1), [1:2), [2:3) [3:inf)]
    if signllr < 50
        binNum = 8;
        if signllr < 3
            binNum = 7;
            if signllr < 2    
                binNum = 6;
                if signllr < 1
                    binNum = 5;
                    if signllr < 0
                        binNum = 4;
                        if signllr < -1
                            binNum = 3;
                            if signllr < -2
                                binNum = 2;
                                if signllr < -3
                                    binNum = 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

