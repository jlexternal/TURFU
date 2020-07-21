% Parameter Finding Script across multiple subjects
%   on the Glaze et al (2015) model + Inference Noise
%
%   Name:   infNoiseGlaze_paramLogger.m
% Author:   Jun Seok Lee
%   Date:   March 2019
%   Type:   Script
%   Note:   Adapted for TURFU experiment

% Excluded subjects:

excluded = [14 20 21 22 27]; % Subject numbers to exclude

% Global variables
global blockFilter; % 'volatility' or 'direction'
global globCond;    % 1 or 2
global globSubject; % subject number reference
global globParticleCount; % amount of particles to use

% Settings
subjects = [28:30];

globParticleCount = 5000;

% Data structure for parameter storage
infNoiseParamsLog = struct;
% The long-ass loop
for i = subjects
    % skip excluded subject numbers
    if ismember(i,excluded)
        infNoiseParamsLog(i).past    = [NaN NaN];
        infNoiseParamsLog(i).future  = [NaN NaN];
        infNoiseParamsLog(i).low     = [NaN NaN];
        infNoiseParamsLog(i).high	 = [NaN NaN];
        continue;
    end
    
    globSubject = i; % specify subject number
    
    % specify the 4 conditions (2 directions x 2 volatilities) and find params
    for filter = {'direction'} % 'volatility'}
        blockFilter = filter;
        for j = 1:2
            globCond = j;
            infNoiseGlaze_paramSearchBADS;
            % Log the parameters in the correct field in the data structure
            if strcmpi(filter,'direction')
                if j == 1
                    infNoiseParamsLog(i).past   = optiParams;
                else
                    infNoiseParamsLog(i).future = optiParams;
                end
            elseif strcmpi(filter,'volatility')
                if j == 1
                    infNoiseParamsLog(i).low = optiParams;
                else
                    infNoiseParamsLog(i).high = optiParams;
                end
            else
                error('Check the values for the global conditions. Something went wrong.');
            end
        end
    end
    save('infNoiseParamsLog.mat','infNoiseParamsLog');
end

