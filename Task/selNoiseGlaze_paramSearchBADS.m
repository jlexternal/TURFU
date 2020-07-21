% Parameter Search with *BADS* on Sequential Particle Filter 
%   on the Glaze et al (2015) model + Inference Noise + Selection Noise
% Author:   Jun Seok Lee
% Date:     February 2019
%
% *BADS: Bayesian Adaptive Direct Search of Acerbi and Ma, 2017

% Global variables
global blockFilter;
global globCond;
global globSubject;

if strcmpi(blockFilter, 'direction') & globCond == 1
    condStr = 'postdictive condition.';
elseif strcmpi(blockFilter, 'direction') & globCond == 2
    condStr = 'predictive condition.';
elseif strcmpi(blockFilter, 'volatility') & globCond == 1
    condStr = 'low volatility condition.';
elseif strcmpi(blockFilter, 'volatility') & globCond == 2
    condStr = 'high volatility condition.';
end

elapsed = tic;
% Load subject file
subject = globSubject; % choose from 1 through #
filename = dir(sprintf('Data/TURFU_S%02d_*.mat', subject));
filename = filename.name;
load(sprintf('Data/%s',filename));

disp(['Analyzing subject ' num2str(subject) ' on blocks of ' condStr])

addpath(genpath('bads-master')); % load the BADS files into working directory

% Given that the objective function is the log-likelihood of REJECTED particles,
%  we do not have to take the negative of said log-likelihood as per the instructions
%  of Acerbi on his github page.

% Search Parameters
%   Given that sigma and hazard rate are our parameters, the vectors below hold
%   values in the order: [hazRate sigma eta]

startPt = [rand()./2 rand() rand()]; % vector corresponding to values of parameters to be estimated
              %     Side note: should probably set this to various functions on the
              %     parameter map to see if the optimum is found around the same
              %     values on the noisy surface
lBound  = [0 0 0]; % HARD lower bound of parameter values
uBound  = [1 2 2]; % HARD upper bound of parameter values
pLBound = [.05 0 0]; % Plausible lower bound of parameter values (optional)
pUBound = [.9 1 1]; % Plausible upper bound of parameter values (optional)

% Changing parameter tolerance options on BADS
options = [];
options.TolMesh = 0.01;

% Initializing values for visuals
txt = 1;
x1 = [];
clf;

% BADS
[optiParams, fVal] = bads(@selNoiseGlaze_PF, startPt, lBound, uBound, pLBound, pUBound, [], options)

toc(elapsed);

