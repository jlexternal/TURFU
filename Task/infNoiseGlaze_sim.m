% Model simulator script for the Glaze (2015) model + Inference Noise
% 
%   Name:   infNoiseGlaze_sim.m
% Author:   Jun Seok Lee
%   Date:   April 2018
%   Type:   Script
%   Note:   Adapted for TURFU experiment

% Description:
% This script can be used to generate instances of the model given specified
% parameters, which it can then also recover if paramRecovery is set to true.

clearvars -except inclRespRate meanRespRate glazeOrigMeanRespRate glazeOrigInclRespRate;

% Global variables
global blockFilter;
global globParticleCount;

% Load parameters obtained from model fitting
load('infNoiseParamsLog','infNoiseParamsLog');

% BADS settings
addpath(genpath('bads-master')); % load the BADS files into working directory

% Search Parameters
% [hazRate sigma]
startPt = [rand() rand()]; % vector corresponding to values of parameters to be estimated
lBound  = [0 0]; % HARD lower bound of parameter values
uBound  = [1 2]; % HARD upper bound of parameter values
pLBound = [.3 0]; % Plausible lower bound of parameter values (optional)
pUBound = [.9 1]; % Plausible upper bound of parameter values (optional)
% Changing parameter tolerance options on BADS
options = [];
options.TolMesh = 0.01;
% /BADS settings

% Static variables (set and forget)
n_subjects  = 30;   % set to number of subjects we are analyzing
excluded    = [21]; % excluded subject numbers

% Dynamic variables (CHANGE VALUES HERE BASED ON ANALYSIS)
blockFilter     = 'direction';  % 'direction' or 'volatility'
nparticles      = 5000;         % number of times to be simulated
globParticleCount = nparticles;
paramRecovery   = false;        % set to true if parameters are to be recovered
exclude         = true;         % set to true to exclude shitty results
specSubjects    = false;        % set to true if analyzing only certain subjects
condition      = [2];          % set to 1 if for 1st condition, 2 for 2nd, 1:2 for both

if exclude
    excluded = [14 20 21 22 27];
end

if specSubjects
    subjects = [1];          % specify subjects here
else
    subjects = [1:n_subjects];
end

% Data structures to hold simulated behavior
simBehavior = struct;

% stuff for model fit to data, specifically p(choice2)|seqllr
infnoisechoice1     = zeros(30,8);
infnoisechoicetot   = zeros(30,8);

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
        disp(['Analyzing subject ' num2str(isubj) ' on condition ' num2str(cond) '...'])
        
        % Reference for condition and fitted parameter
        if strcmpi(blockFilter, 'direction')
            blocks = find([expe.blck.taskid] == cond);	% taskid refers to inference direction in exp.
            haz     = infNoiseParamsLog(isubj).past(1);
            sigm    = infNoiseParamsLog(isubj).past(2);
        elseif strcmpi(blockFilter, 'volatility')
            blocks = find([expe.blck.condtn] == cond);	% condtn refers to volatility of block
            if cond == 1
                haz     = infNoiseParamsLog(isubj).low(1);
                sigm    = infNoiseParamsLog(isubj).low(2);
            else
                haz     = infNoiseParamsLog(isubj).high(1);
                sigm    = infNoiseParamsLog(isubj).high(2);
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
            for iparticle = 1:nparticles
                if strcmpi(blockFilter, 'direction') && condition == 1
                [simBehavior(isubj).rslt(block).resp(iparticle,:), simBehavior(isubj).rslt(block).belief(iparticle,:)] = ...
                                                                    infNoiseGlaze_simulator(haz, sigm, seqllrs, nsamples,0);
                elseif strcmpi(blockFilter, 'direction') && condition == 2
                    [simBehavior(isubj).rslt(block).resp(iparticle,:), simBehavior(isubj).rslt(block).belief(iparticle,:)] = ...
                                                                    infNoiseGlaze_simulator(haz, sigm, seqllrs, nsamples,1);
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
                tempChoicetot = infnoisechoicetot(isubj,binIndex(expe.blck(block).seqllr(itrial)));
                infnoisechoicetot(isubj,binIndex(expe.blck(block).seqllr(itrial))) = tempChoicetot + nparticles;
                
                tempChoice1 = infnoisechoice1(isubj,binIndex(expe.blck(block).seqllr(itrial))); % get previous value in bin
                tempTrialChoice1 = simBehavior(isubj).rslt(block).resp(:,itrial+1);     % extract array from structure for easiness
                infnoisechoice1(isubj,binIndex(expe.blck(block).seqllr(itrial))) = tempChoice1 + numel(tempTrialChoice1(tempTrialChoice1==1)); % update bin
            end
            
            % model log-likelihood calculation
            corTrls = simBehavior(isubj).rslt(block).resp(2:73)-expe.rslt(block).resp(2:73); % correct generations
            modelL(isubj,iblock) = numel(corTrls(corTrls==0))/72; % take of likelihood of model to behavior on block
        
            iblock = iblock + 1;
        end
        
        if paramRecovery
            % Initializing values for BADS visuals
            txt = 1;
            x1 = [];
            clf;
            % run BADS
            [optiParams, fVal] = bads(@infNoiseGlaze_PF, startPt, lBound, uBound, pLBound, pUBound, [], options)

            % Register the recovered parameters into the original paramsLog structure
            if strcmpi(blockFilter, 'direction')
                if cond == 1
                    infNoiseParamsLog(isubj).pastRecd   = optiParams;
                else
                    infNoiseParamsLog(isubj).futureRecd = optiParams;
                end
            elseif strcmpi(blockFilter, 'volatility')
                if cond == 1
                    infNoiseParamsLog(isubj).lowRecd    = optiParams;
                else
                    infNoiseParamsLog(isubj).highRecd   = optiParams;
                end
            end
        end
    end
end

% save simulated choices and beliefs to file
dateStr = datestr(datenum(date),'ddmmyyyy');
if strcmpi(blockFilter,'direction')
    if condition == 1
        save(['infNoiseGlaze_sim_post_struct_' dateStr], 'simBehavior');
    elseif condition == 2
        save(['infNoiseGlaze_sim_pred_struct_' dateStr], 'simBehavior');
    end
end

% Calculate total likelihood that the data came from the model
totalL = [];
for isubj = 1%:numel(modelL(:,1))
    if modelL(isubj,1) == 0
        continue;
    end
    % list the mean likelihood across the blocks for each subject
    totalL = vertcat(totalL, sum(modelL(isubj,:)));%/numel(modelL(isubj,:))); 
end
totalL  = mean(totalL); % take mean likelihood across the subjects
totalLL = log(totalL);  % take the log of the mean likelihood

%% Convert to percentages and plot

infNoiseRespRate = infnoisechoice1./infnoisechoicetot;
exclRows = isnan(infNoiseRespRate(:,1));
infNoiseInclRespRate = infNoiseRespRate(~exclRows,:);
infNoiseMeanRespRate = mean(infNoiseInclRespRate);

scatter([1:8], infNoiseMeanRespRate);
hold on;
errorbar([1:8], infNoiseMeanRespRate, std(infNoiseInclRespRate));
hold off;

figure(2);
hold on;
for itest = 1:25
    scatter([1:8],infNoiseInclRespRate(itest,:));
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

