% Parameter Search with *BADS* on the Original Glaze (2015) model
%
%   Name:   glazeOriginal_paramSearchBADS
%   Type:   Script
% Author:   Jun Seok Lee
%   Date:   March 2019
%   Note:   Adapted for TURFU
%
% *BADS: Bayesian Adaptive Direct Search of Acerbi and Ma, 2017
%
%   NOTE: This code is modified for the data structures created by the TURFU experiment
% Parameter Search with *BADS* on the Determinstic Glaze model
%   on the Glaze et al (2015) model
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

% Search Parameters
%   Given that  hazard rate is our parameter, the vectors below hold
%   values in the order: [hazRate]

startPt = [rand().*.5]; % vector corresponding to values of parameters to be estimated
              %     Side note: should probably set this to various functions on the
              %     parameter map to see if the optimum is found around the same
              %     values on the noisy surface
lBound  = [0]; % HARD lower bound of parameter values
uBound  = [.51]; % HARD upper bound of parameter values
pLBound = [.01]; % Plausible lower bound of parameter values (optional)
pUBound = [.5]; % Plausible upper bound of parameter values (optional)


% Changing parameter tolerance options on BADS
options = [];
options.TolMesh = 0.01;

% Initializing values for visuals
txt = 1;
x1 = [];
clf;

% BADS
[optiParams, fVal] = bads(@glazeOriginal, startPt, lBound, uBound, pLBound, pUBound, [], options)


toc(elapsed);


