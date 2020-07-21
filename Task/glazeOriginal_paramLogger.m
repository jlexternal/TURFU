% Parameter Finding Script across multiple subjects
%   on the Glaze et al (2015) model
%
%   Name:   glazeOriginal_paramLogger.m
% Author:   Jun Seok Lee
%   Date:   April 2019
%   Type:   Script
%   Note:   Adapted for TURFU experiment


% Global variables
global blockFilter; % 'volatility' or 'direction'
global globCond;    % 1 or 2
global globSubject; % subject number reference
global globParticleCount; % amount of particles to use

% Static variables (set and forget)
n_subjects  = 30;   % set to number of subjects we are analyzing
excluded    = [14 20 21 22 27]; % excluded subject numbers

globParticleCount = 2000;

% Data structure for parameter storage
glazeOrigParamLog = struct;

% Parameter search and log
for i = 1:n_subjects
    % skip excluded subject numbers
    if ismember(i,excluded)
        glazeOrigParamsLog(i).past    = [NaN];
        glazeOrigParamsLog(i).future  = [NaN];
        glazeOrigParamsLog(i).low     = [NaN];
        glazeOrigParamsLog(i).high	  = [NaN];
        continue;
    end
    
    globSubject = i; % specify subject number
    
    % specify the 4 conditions (2 directions x 2 volatilities) and find params
    for filter = {'direction'} %'volatility'}
        blockFilter = filter;
        for j = 1:2
            globCond = j;
            glazeOriginal_paramSearchBADS;
            % Log the parameters in the correct field in the data structure
            if strcmpi(filter,'direction')
                if j == 1
                    glazeOrigParamsLog(i).past  = optiParams;
                else
                    glazeOrigParamsLog(i).future= optiParams;
                end
            elseif strcmpi(filter,'volatility')
                if j == 1
                    glazeOrigParamsLog(i).low = optiParams;
                else
                    glazeOrigParamsLog(i).high = optiParams;
                end
            else
                error('Check the values for the global conditions. Something went wrong.');
            end
        end
    end
    %save('glazeOrigParamLog.mat', 'glazeOrigParamLog');
end

