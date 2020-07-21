% Parameter Search with *BADS* on Sequential Particle Filter 
%   on the Forward Glaze model
% Author:   Jun Seok Lee
% Date:     February 2019
%
% *BADS: Bayesian Adaptive Direct Search of Acerbi and Ma, 2017

clf;
% Load subject file
%   global globSubj;
%   subject = globSubj;
subject = 1;

% Temporary code for testing:
global blockFilter;
global globCond;

blockFilter = 'condition';
globCond = 2;

elapsed = tic;
filename = sprintf('S%02d/expe.mat', subject);
load(filename);

addpath(genpath('bads-master')); % load the BADS files into working directory

% Given that the objective function is the log-likelihood of REJECTED particles,
%  we do not have to take the negative of said log-likelihood as per the instructions
%  of Acerbi on his github page.

% Search Parameters
%   Given that sigma and hazard rate are our parameters, the vectors below hold
%   values in the order: [sigma hazRate]

startPt = [rand() rand()]; % vector corresponding to values of parameters to be estimated
              %     Side note: should probably set this to various functions on the
              %     parameter map to see if the optimum is found around the same
              %     values on the noisy surface
lBound  = [0 0]; % HARD lower bound of parameter values
uBound  = [1 1]; % HARD upper bound of parameter values
pLBound = [.1 .01]; % Plausible lower bound of parameter values (optional)
pUBound = [.9 .5]; % Plausible upper bound of parameter values (optional)

% To use parallel processing at the particle filter level, set to true
parProc = false;

% Changing parameter tolerance options on BADS
options = [];
options.TolMesh = 0.01;

% Initializing values for visuals
txt = 1;
x1 = [];

% BADS
if ~parProc
    [optiParams, fVal] = bads(@glazeForward_PF, startPt, lBound, uBound, pLBound, pUBound, [], options)
else
    %[optiParams, fVal] = bads(@seqPartFilt, startPt, lBound, uBound, pLBound, pUBound, [], options)
end

toc(elapsed);

