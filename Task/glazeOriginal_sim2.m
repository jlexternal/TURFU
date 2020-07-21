% Model simulator script for the Glaze (2015) model 
%
%   Name:   glazeOriginal_sim2.m
% Author:   Jun Seok Lee
%   Date:   April 2018
%   Type:   Script
%   Note:   Adapted for the TURFU experiment.
%           This script is different from glazeOriginal_sim.m in that it is not
%           recovering fitted parameters on subjects, but rather, from a set of
%           predetermined parameters, it generates behavior based on some chosen
%           experimental sequences from subjects.


% Load parameters obtained from model fitting
load('glazeOrigParamsLog','glazeOrigParamsLog');
clearvars -except glazeOrigParamsLog;

%%%%%%%%%%%%%%%%%%%%%%%%% BADS settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('bads-master')); % load the BADS files into working directory
% Search Parameters
%   Given that sigma and hazard rate are our parameters, the vectors below hold
%   values in the order: [hazRate]
startPt = [rand().*.5]; % vector corresponding to values of parameters to be estimated
              %     Side note: should probably set this to various functions on the
              %     parameter map to see if the optimum is found around the same
              %     values on the noisy surface
lBound  = [0]; % HARD lower bound of parameter values
uBound  = [.7]; % HARD upper bound of parameter values
pLBound = [.01]; % Plausible lower bound of parameter values (optional)
pUBound = [.69]; % Plausible upper bound of parameter values (optional)
% Changing parameter tolerance options on BADS
options = [];
options.TolMesh = 0.01;
%%%%%%%%%%%%%%%%%%%%%%%% /BADS settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Static variables (set and forget)
paramSet    = [.05:.05:.6];	% H parameters to simulate
expeSet     = [1:8];       	% subject numbers

%% Construct experimental structure and simulated behavior

% create simBehavior structure with null data
simBehavior = struct;
for i = 1:numel(paramSet)
    simBehavior(i).H_val        = paramSet(i);
    simBehavior(i).expe         = struct;
    simBehavior(i).rslt         = struct;
    simBehavior(i).H_val_recd   = zeros(1,numel(expeSet));
    for j = 1:numel(expeSet)
        % fill in the experiment sequence data
        simBehavior(i).expe(j).seqllr   = zeros(1,72);
        % create null structure to be filled in w/ responses later on
        simBehavior(i).rslt(j).resp = zeros(1,73);
    end
end

% fill in the expe structure for each H value with experimental llrs
for expeIndex = expeSet
    filename = dir(sprintf('Data/TURFU_S%02d_*.mat', expeIndex));
    filename = filename.name;
    load(sprintf('Data/%s',filename));
    
    for hIndex = 1:numel(paramSet)
        % Note: expIndex+2 is used to get the 8 real blocks
        simBehavior(hIndex).expe(expeIndex).seqllr = expe.blck(expeIndex+2).seqllr;
    end
end

% simulate behavior based on setup above
for hIndex = 1:numel(paramSet)
    disp(num2str(paramSet(hIndex)));
    for expeIndex = 1:numel(expeSet)
        seqllrs = simBehavior(hIndex).expe(expeIndex).seqllr; 
        hazR = paramSet(hIndex);
        simBehavior(hIndex).rslt(expeIndex).resp = glazeOriginal_simulator(hazR, seqllrs);
        
        % Initializing values for BADS visuals
        txt = 1;
        x1 = [];
        clf;
        % run BADS
        [optiParams, fVal] = bads(@glazeOriginal_v_modelRec, startPt, lBound, uBound, pLBound, pUBound, [], options);
        
        simBehavior(hIndex).H_val_recd(expeIndex) = optiParams;
    end
end
