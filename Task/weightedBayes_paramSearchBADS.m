% Parameter Search with *BADS* on the Weighted Bayes LLR model
%
% Type:     Script
% Author:   Jun Seok Lee
% Date:     February 2019
%
% *BADS: Bayesian Adaptive Direct Search of Acerbi and Ma, 2017
clf;

global globSubj;

elapsed = tic;
% Load subject file
subject = globSubj;

if modelRecovery == false
    filename = sprintf('S%02d/expe.mat', subject);
else
    filename = sprintf('S%02d_m/expe.mat', subject);
end
load(filename);

addpath(genpath('bads-master')); % load the BADS files into working directory

% Search Parameters
%   Given that sigma and hazard rate are our parameters, the vectors below hold
%   values in the order: [sigma hazRate]

startPt = [rand()]; % vector corresponding to values of parameters to be estimated
              %     Side note: should probably set this to various functions on the
              %     parameter map to see if the optimum is found around the same
              %     values on the noisy surface
lBound  = [0]; % HARD lower bound of parameter values
uBound  = [1]; % HARD upper bound of parameter values
pLBound = [.01]; % Plausible lower bound of parameter values (optional)
pUBound = [.99]; % Plausible upper bound of parameter values (optional)


% Initializing values for visuals
txt = 1;
x1 = [];

% BADS
[optiParams, fVal] = bads(@weightedBayes, startPt, lBound, uBound, pLBound, pUBound);


toc(elapsed);

