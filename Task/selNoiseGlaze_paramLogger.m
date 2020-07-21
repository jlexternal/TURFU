% Parameter Finding Script across multiple subjects
%   on the Glaze et al (2015) model + Inference Noise + Selection Noise
%
%   Name:   selNoiseGlaze_paramLogger.m
% Author:   Jun Seok Lee
%   Date:   April 2019
%   Type:   Script
%   Note:   Adapted for TURFU experiment
%           *There are three parameters to search amongst: H, sigma, eta

% Excluded subjects:

excluded = [14 20 21 22 27]; % Subject numbers to exclude

% Global variables
global blockFilter; % 'volatility' or 'direction'
global globCond;    % 1 or 2
global globSubject; % subject number reference
global globParticleCount; % amount of particles to use

globParticleCount = 2000;

% Data structure for parameter storage
%selNoiseParamsLog = struct;
% The long-ass loop
for i = 1:30
    % skip excluded subject numbers
    if ismember(i,excluded)
        selNoiseParamsLog(i).past    = [NaN NaN NaN];
        selNoiseParamsLog(i).future  = [NaN NaN NaN];
        selNoiseParamsLog(i).low     = [NaN NaN NaN];
        selNoiseParamsLog(i).high	 = [NaN NaN NaN];
        continue;
    end
    
    globSubject = i; % specify subject number
    
    % specify the 4 conditions (2 directions x 2 volatilities) and find params
    for filter = {'direction'}% 'volatility'}
        blockFilter = filter;
        for j = 2
            globCond = j;
            selNoiseGlaze_paramSearchBADS;
            % Log the parameters in the correct field in the data structure
            if strcmpi(filter,'direction')
                if j == 1
                    selNoiseParamsLog(i).past   = optiParams;
                else
                    selNoiseParamsLog(i).future = optiParams;
                end
            elseif strcmpi(filter,'volatility')
                if j == 1
                    selNoiseParamsLog(i).low    = optiParams;
                else
                    selNoiseParamsLog(i).high   = optiParams;
                end
            else
                error('Check the values for the global conditions. Something went wrong.');
            end
        end
    end
    save('selNoiseParamsLog.mat','selNoiseParamsLog');
end

