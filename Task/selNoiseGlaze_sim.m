% Model simulator script for the Glaze (2015) model + Inference Noise + Selection
% Noise

%   Name:   selNoiseGlaze_sim.m
% Author:   Jun Seok Lee
%   Date:   April 2018
%   Type:   Script
%   Note:   Adapted for TURFU experiment

% Description: 
% This script can be used to generate instances of the model given specified
% parameters, which it can then also recover if paramRecovery is set to true.

clearvars -except inclRespRate meanRespRate glazeOrigMeanRespRate glazeOrigInclRespRate infNoiseMeanRespRate infNoiseInclRespRate;

% Global variables
global blockFilter;
global globParticleCount;

% Load parameters obtained from model fitting
load('selNoiseParamsLog','selNoiseParamsLog');

% BADS settings
addpath(genpath('bads-master')); % load the BADS files into working directory

% Search Parameters
% [hazRate sigma eta]
startPt = [rand().*.5 rand() rand()]; % vector corresponding to values of parameters to be estimated
lBound  = [0   0  0]; % HARD lower bound of parameter values
uBound  = [1   1  2]; % HARD upper bound of parameter values
pLBound = [.05 0  .01]; % Plausible lower bound of parameter values (optional)
pUBound = [.9  .9 1]; % Plausible upper bound of parameter values (optional)
% Changing parameter tolerance options on BADS
options = [];
options.TolMesh = 0.01;
% /BADS settings

% Static variables (set and forget)
n_subjects  = 30;   % set to number of subjects we are analyzing
excluded    = [21]; % excluded subject numbers

if n_subjects ~= numel(selNoiseParamsLog) % error catching
    error('Number of subjects is not equal to the number of entries in the parameters log.');
end

% Dynamic variables (CHANGE VALUES HERE BASED ON ANALYSIS)
blockFilter     = 'direction';  % 'direction' or 'volatility'
nparticles      = 5000;          % number of times to be simulated
globParticleCount = nparticles;
paramRecovery   = false;        % set to true if parameters are to be recovered
exclude         = true;         % set to true to exclude shitty results
specSubjects    = false;        % set to true if analyzing only certain subjects
condition      = [2];          % set to 1 if for 1st condition, 2 for 2nd : RUN ONLY ONE AT A TIME

if exclude
    excluded = [14 20 21 22 27];
end
if specSubjects
    subjects = [1 2];          % specify subjects here
else
    subjects = [1:n_subjects];
end

% Data structures to hold simulated behavior
simBehavior = struct;

% stuff for model fit to data, specifically p(choice2)|seqllr
selnoisechoice1     = zeros(30,8);
selnoisechoicetot   = zeros(30,8);

modelLL = zeros(30); % structure to hold model likelihoods to subject data

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
    disp(['running subject ' num2str(isubj) '...']);
    
    % identify condition that parameters were found in
    for cond = condition
        
        % Block and parameter choice (based on analysis)
        if strcmpi(blockFilter, 'direction')
            blocks = find([expe.blck.taskid] == cond);	% taskid refers to inference direction in exp.
            haz     = selNoiseParamsLog(isubj).past(1);
            sigm    = selNoiseParamsLog(isubj).past(2);
            eta     = selNoiseParamsLog(isubj).past(3);
        elseif strcmpi(blockFilter, 'volatility')
            blocks = find([expe.blck.condtn] == cond);	% condtn refers to volatility of block
            if cond == 1
                haz     = selNoiseParamsLog(isubj).low(1);
                sigm    = selNoiseParamsLog(isubj).low(2);
                eta     = selNoiseParamsLog(isubj).low(3);
            else
                haz     = selNoiseParamsLog(isubj).high(1);
                sigm    = selNoiseParamsLog(isubj).high(2);
                eta     = selNoiseParamsLog(isubj).low(3);
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
            
            % simulate behavior based on params for block and store into local memory
            for iparticle = 1:nparticles
                if strcmpi(blockFilter, 'direction') && condition == 1
                    [simBehavior(isubj).rslt(block).resp(iparticle,:), simBehavior(isubj).rslt(block).belief(iparticle,:)] = ...
                                                                            selNoiseGlaze_simulator(haz, sigm, eta, seqllrs, nsamples, 0);
                elseif strcmpi(blockFilter, 'direction') && condition == 2
                    [simBehavior(isubj).rslt(block).resp(iparticle,:), simBehavior(isubj).rslt(block).belief(iparticle,:)] = ...
                                                                            selNoiseGlaze_simulator(haz, sigm, eta, seqllrs, nsamples, 1);
                end
            end
            
            if paramRecovery
                % RISKY CODE (need to verify that this does not change the save file)
                % replace human responses with simulated responses directly into expe variable
                expe.rslt(block).resp = simBehavior(isubj).rslt(block).resp; 
            end
            
            
            % testing model fit to data stuff
            for itrial = 1:72
                
                % update number of total particles for the seqllr bin on the trial
                tempChoicetot = selnoisechoicetot(isubj,binIndex(expe.blck(block).seqllr(itrial)));
                selnoisechoicetot(isubj,binIndex(expe.blck(block).seqllr(itrial))) = tempChoicetot + nparticles;
                
                tempChoice1 = selnoisechoice1(isubj,binIndex(expe.blck(block).seqllr(itrial))); % get previous value in bin
                tempTrialChoice1 = simBehavior(isubj).rslt(block).resp(:,itrial+1);     % extract array from structure for easiness
                selnoisechoice1(isubj,binIndex(expe.blck(block).seqllr(itrial))) = tempChoice1 + numel(tempTrialChoice1(tempTrialChoice1==1)); % update bin
            end
            
            % model log-likelihood calculation
            corTrls = simBehavior(isubj).rslt(block).resp(2:73)-expe.rslt(block).resp(2:73); % correct generations
            modelL(isubj,iblock) = numel(corTrls(corTrls==0))/72; % take likelihood of model to behavior on block
        
            iblock = iblock + 1;
        end
        
        if paramRecovery
            % Initializing values for BADS visuals
            txt = 1;
            x1 = [];
            clf;
            % run BADS
            [optiParams, fVal] = bads(@selNoiseGlaze_PF, startPt, lBound, uBound, pLBound, pUBound, [], options)

            % Register the recovered parameters into the original selNoiseParamsLog structure
            if strcmpi(blockFilter, 'direction')
                if cond == 1
                    selNoiseParamsLog(isubj).pastRecd   = optiParams;
                else
                    selNoiseParamsLog(isubj).futureRecd = optiParams;
                end
            elseif strcmpi(blockFilter, 'volatility')
                if cond == 1
                    selNoiseParamsLog(isubj).lowRecd    = optiParams;
                else
                    selNoiseParamsLog(isubj).highRecd   = optiParams;
                end
            end
        end
    end
end

% save simulated choices and beliefs to file
dateStr = datestr(datenum(date),'ddmmyyyy');

if strcmpi(blockFilter,'direction')
    if condition == 1
        save(['selNoiseGlaze_sim_post_struct_' dateStr], 'simBehavior');
    elseif condition == 2
        save(['selNoiseGlaze_sim_pred_struct_' dateStr], 'simBehavior');
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

selNoiseRespRate = selnoisechoice1./selnoisechoicetot;
exclRows = isnan(selNoiseRespRate(:,1));
selNoiseInclRespRate = selNoiseRespRate(~exclRows,:);
selNoiseMeanRespRate = mean(selNoiseInclRespRate);

scatter([1:8], selNoiseMeanRespRate);
hold on;
errorbar([1:8], selNoiseMeanRespRate, std(selNoiseInclRespRate));
hold off;

figure(2);
hold on;
for itest = 1:25
    scatter([1:8],selNoiseInclRespRate(itest,:));
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


