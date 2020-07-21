% Correction on subject mis-mapping blocks

% Run this script to create data sets that are corrected for when subjects switched
% the mapping on certain blocks of the experiment. 

%   Name:   mismapCorrection.m
%   Type:   Script
% Author:   Jun Seok Lee
%   Date:   March 2019

% Subjects and their mismapped blocks (blocks referred to absolutely)


% Statics
n_trials = 73; % switch the useless 1st value
datapath = ['Data'];

% Data structures

subjs = [10]; % specify in the array the subject numbers of participants who made a mistake

blocks = {[8]}; % specify in the structure the block numbers that need color mapping switches

isubj = 1;
for subj = subjs
    
    %load subject information
    subject = subj;
    filename = dir(sprintf('Data/TURFU_S%02d_*.mat', subject));
    filename = filename.name;
    load(sprintf('Data/%s',filename));
    
    jblck = 1;
    for j = blocks{isubj} % j is the block number
            for k = 1:n_trials
                if expe.rslt(blocks{isubj}(jblck)).resp(k) == 1
                    expe.rslt(blocks{isubj}(jblck)).resp(k) = 2;
                else
                    expe.rslt(blocks{isubj}(jblck)).resp(k) = 1;
                end
            end
        jblck = jblck + 1;
    end
    
    %output fixed .mat file
    filename = filename(1:end-4);
    filename = strcat(filename, '_corrected.mat');
    filename = fullfile(datapath,filename);
    save(filename,'expe');
    isubj = isubj + 1;
end
