% Model-based analysis on TURFU subjects and models
%
%   Name:   modelBased_org.m
%   Type:   Script
% Author:   Jun Seok Lee
%   Date:   May 2019
% 

% In order for this code to function and do a proper analysis, the simulations should
% have already been run and stored into a saved file. 

clf;
close all;

% colors for plot
postRGB = [1 .75 .5];
predRGB = [.75 .75 1];

%% Setup
% Static variables (set and forget)
n_subjects  = 30;   % set to number of subjects we are analyzing
excluded    = [21]; % excluded subject numbers

% **Dynamic variables (CHANGE VALUES HERE BASED ON ANALYSIS)**
blockFilter     = 'direction'; % 'direction' or 'volatility'
condition       = 1;           % 1 for past/low; 2 for future/high (direction/volatility)
befAftTrials    = 4;           % how many trials before and after change point
exclude         = true;        % set to true to exclude shitty results
if exclude
    excluded = [14 20 21 22 27];
end
% set below values to 'false' if skipping a certain model for whatever reason
model_det   = true;     % deterministic model
model_inf   = true;     % inference noise model
model_sel   = true;     % inf. noise + selection noise model
% **                                                        **

if condition == 1
    condStr = 'post';
elseif condition == 2
    condStr = 'pred';
end

% Plot legend text
if strcmpi(blockFilter, 'direction')
    txtCond1 = 'Postdiction <-';
    txtCond2 = 'Prediction  ->';
else
    txtCond1 = 'Low volatility';
    txtCond2 = 'High volatility';
end

% Model statistics structures 
load('modelStruct','modelStruct');

% Note: confidence cutoff thresholds are stored in confCutoff.mat in format
% (condition, subject number, 1:negllr/2:posllr)

%% Data filtering and organization

for isubj = 1:n_subjects
    %disregard excluded subjects
    if ismember(isubj, excluded)
       continue; 
    end
    
    %load subject information
    subject  = isubj;
    filename = dir(sprintf('Data/TURFU_S%02d_*.mat', subject));
    filename = filename.name;
    load(sprintf('Data/%s',filename));
    
    disp(['Analyzing subject ' num2str(isubj) ' on condition ' num2str(condition)]);
    
    % Analyze deterministic model simulation results
    if model_det
        disp('on deterministic model');
        %load saved simStructure for this model
        filename = dir(sprintf('SavedStructures/glazeOriginal_sim_%s_struct_*.mat',condStr));
        filename = filename.name;
        load(sprintf('SavedStructures/%s',filename));
        
        tempRevRate = [];
        tempRevConf = [];
        
        temp_binReps = zeros(1,befAftTrials.*2);
        temp_binTots = zeros(1,befAftTrials.*2);
        temp_signConf = zeros(1,befAftTrials.*2);
        
        %go through the blocks
        for iblock = 3:10 
            % analyze blocks corresponding to the specific block condition
            if expe.blck(iblock).taskid ~= condition
                continue; 
            end
            
            %go through the trials (itrial is the true experimental trial number)
            for itrial = 2:72 
                
            % Reversal Curve analysis
                if expe.blck(iblock).seqdir(itrial) ~= expe.blck(iblock).seqdir(itrial-1) && itrial~=2  % identify trial after a reversal
                    reversalRange = itrial-befAftTrials+1:itrial+befAftTrials;            % trial numbers to consider and record
                    tempRevRate = vertcat(tempRevRate, ...
                                          simBehavior(isubj).rslt(iblock).resp(reversalRange)==expe.blck(iblock).seqdir(itrial));
                                      
                    tempRevConf = vertcat(tempRevConf, simBehavior(isubj).rslt(iblock).conf(reversalRange-1)-1);
                end
                
            % Repetition Rate on Signed (by prev choice) Evidence analysis
                simChoices = (simBehavior(isubj).rslt(iblock).resp(itrial)-1.5)*-2;   % convert prev choice to sign (+/- 1)
                tempRepeats = simBehavior(isubj).rslt(iblock).resp(itrial+1) == ...   % find repeated choice trial
                              simBehavior(isubj).rslt(iblock).resp(itrial);
                tempRepeats = double(tempRepeats);                                    % convert logical -> numerical
                iseqllr = expe.blck(iblock).seqllr(itrial);                           % get sequence llr
                temp_binIndex = binIndex(iseqllr.*simChoices);                        % identify proper total signed llr bin
                temp_binTots(temp_binIndex) = temp_binTots(temp_binIndex)+1;          % add to total signed llr encountered
                if tempRepeats == 1                                                   % if particle repeated, count to proper llr bin
                    temp_binReps(temp_binIndex) = temp_binReps(temp_binIndex) + 1;    % add to the repeat bin
                end
                if simBehavior(isubj).rslt(iblock).conf(itrial) == 2                    % if high confidence
                    temp_signConf(temp_binIndex) = temp_signConf(temp_binIndex) + 1;      % add to appopriate bin
                end
                
            end %loop over trials
        end %loop over blocks
        
        % update holding structures
            % 1 reversal analysis
            % 2 repeat/stay analysis
        tempRepRate = temp_binReps./temp_binTots;
        tempRepRate(isnan(tempRepRate)) = 0;
        if condition == 1
            modelStruct(isubj).det.past.rev      = sum(tempRevRate)./size(tempRevRate,1); %1
            modelStruct(isubj).det.past.revconf  = sum(tempRevConf)./size(tempRevConf,1); %1
            modelStruct(isubj).det.past.rep = tempRepRate; %2
            modelStruct(isubj).det.past.signconf = temp_signConf./temp_binTots; %2 high confidence rate on signed llr
        elseif condition == 2
            modelStruct(isubj).det.future.rev       = sum(tempRevRate)./size(tempRevRate,1); %1
            modelStruct(isubj).det.future.revconf   = sum(tempRevConf)./size(tempRevConf,1); %1
            modelStruct(isubj).det.future.rep       = tempRepRate; %2
            modelStruct(isubj).det.future.signconf  = temp_signConf./temp_binTots; %2 high confidence rate on signed llr
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Analyze inference noise model simulation results
    if model_inf
        disp('on inf noise model');
        %load saved simStructure for this model
        filename = dir(sprintf('SavedStructures/infNoiseGlaze_sim_%s_struct_*.mat',condStr));
        filename = filename.name;
        load(sprintf('SavedStructures/%s',filename));
        
        firstIter = true;
        
        tempRevRate = [];
        tempRevConf = [];
        
        temp_binReps = zeros(1,befAftTrials.*2);
        temp_binTots = zeros(1,befAftTrials.*2);
        temp_signConf = zeros(1,befAftTrials.*2);
        
        %go through the blocks
        for iblock = 3:10 
            % analyze blocks corresponding to the specific block condition
            if expe.blck(iblock).taskid ~= condition
                continue; 
            end
            
            if firstIter % access the number of particles ONCE per subject
                n_particles = numel(simBehavior(isubj).rslt(iblock).resp(:,1));
                firstIter = false;
            end
            
            %%go through the trials (itrial is the true experimental trial number)
            for itrial = 2:72
                
            % Reversal Curve analysis
                if expe.blck(iblock).seqdir(itrial) ~= expe.blck(iblock).seqdir(itrial-1) && itrial~=2   % identify trial after a reversal
                    reversalRange = itrial-befAftTrials+1:itrial+befAftTrials;              % trial numbers to consider and record
                    tempRevRate = vertcat(tempRevRate, ...
                                          simBehavior(isubj).rslt(iblock).resp(:,reversalRange)==expe.blck(iblock).seqdir(itrial));
                    
                    tempRevConf = vertcat(tempRevConf, simBehavior(isubj).rslt(iblock).conf(:,reversalRange-1)-1);
                end
                
            % Repetition Rate on Signed Evidence analysis
                simChoices = (simBehavior(isubj).rslt(iblock).resp(:,itrial)-1.5)*-2;   % convert previous sim choice to sign (+/- 1)
                tempRepeats = simBehavior(isubj).rslt(iblock).resp(:,itrial+1) == ...   % find which particles repeated 
                              simBehavior(isubj).rslt(iblock).resp(:,itrial);
                tempRepeats = double(tempRepeats);                                      % convert logical -> numerical
                iseqllr = expe.blck(iblock).seqllr(itrial);                             % get sequence llr
                for iparticle = 1:n_particles
                    temp_binIndex = binIndex(iseqllr.*simChoices(iparticle));           % identify proper total signed llr bin
                    temp_binTots(temp_binIndex) = temp_binTots(temp_binIndex)+1;        % add to total signed llr encountered
                    if tempRepeats(iparticle) == 1                                      % if particle repeated, count to proper llr bin
                        temp_binReps(temp_binIndex) = temp_binReps(temp_binIndex) + 1;  % add to the repeat bin
                    end
                    if simBehavior(isubj).rslt(iblock).conf(iparticle,itrial) == 2      % if high confidence
                        temp_signConf(temp_binIndex) = temp_signConf(temp_binIndex) + 1;  % add to appopriate bin
                    end
                end
            end %loop over trials
        end %loop over blocks

        % update holding structures
        if condition == 1
            modelStruct(isubj).inf.past.rev         = sum(tempRevRate)./size(tempRevRate,1); %1
            modelStruct(isubj).inf.past.revconf     = sum(tempRevConf)./size(tempRevConf,1); %1
            modelStruct(isubj).inf.past.rep         = temp_binReps./temp_binTots; %2
            modelStruct(isubj).inf.past.signconf    = temp_signConf./temp_binTots; %2 high confidence rate on signed llr
        elseif condition == 2
            modelStruct(isubj).inf.future.rev       = sum(tempRevRate)./size(tempRevRate,1); %1
            modelStruct(isubj).inf.future.revconf   = sum(tempRevConf)./size(tempRevConf,1); %1
            modelStruct(isubj).inf.future.rep       = temp_binReps./temp_binTots; %2
            modelStruct(isubj).inf.future.signconf  = temp_signConf./temp_binTots; %2 high confidence rate on signed llr
        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Analyze inf. noise + selection noise model simulation results
    if model_sel
        disp('on sel noise model');
        %load saved simStructure for this model
        filename = dir(sprintf('SavedStructures/selNoiseGlaze_sim_%s_struct_*.mat',condStr));
        filename = filename.name;
        load(sprintf('SavedStructures/%s',filename));
        
        firstIter = true;
        
        tempRevRate = [];
        tempRevConf = [];
        
        temp_binReps = zeros(1,befAftTrials.*2);
        temp_binTots = zeros(1,befAftTrials.*2);
        temp_signConf = zeros(1,befAftTrials.*2);
        
        %go through the blocks
        for iblock = 3:10 
            % analyze blocks corresponding to the specific block condition
            if expe.blck(iblock).taskid ~= condition
                continue; 
            end
            
            if firstIter % access the number of particles ONCE per subject
                n_particles = numel(simBehavior(isubj).rslt(iblock).resp(:,1));
                firstIter = false;
            end
            
            %go through the trials
            for itrial = 2:72
                
            % Reversal Curve analysis
                if expe.blck(iblock).seqdir(itrial) ~= expe.blck(iblock).seqdir(itrial-1) && itrial~=2    % identify trial after a reversal
                    reversalRange = itrial-befAftTrials+1:itrial+befAftTrials;            % trial numbers to consider and record
                    tempRevRate = vertcat(tempRevRate, ...
                                          simBehavior(isubj).rslt(iblock).resp(:,reversalRange)==expe.blck(iblock).seqdir(itrial));
                                      
                    tempRevConf = vertcat(tempRevConf, simBehavior(isubj).rslt(iblock).conf(:,reversalRange-1)-1);
                end
                
            % Repetition Rate on Signed Evidence analysis
                simChoices = (simBehavior(isubj).rslt(iblock).resp(:,itrial)-1.5)*-2;   % convert previous sim choice to sign (+/- 1)
                tempRepeats = simBehavior(isubj).rslt(iblock).resp(:,itrial+1) == ...   % find which particles repeated 
                              simBehavior(isubj).rslt(iblock).resp(:,itrial);
                tempRepeats = double(tempRepeats);                                      % convert logical -> numerical
                iseqllr = expe.blck(iblock).seqllr(itrial);                             % get sequence llr
                for iparticle = 1:n_particles
                    temp_binIndex = binIndex(iseqllr.*simChoices(iparticle));           % identify proper total signed llr bin
                    temp_binTots(temp_binIndex) = temp_binTots(temp_binIndex)+1;        % add to total signed llr encountered
                    if tempRepeats(iparticle) == 1                                      % if particle repeated, count to proper llr bin
                        temp_binReps(temp_binIndex) = temp_binReps(temp_binIndex) + 1;  % add to the repeat bin
                    end
                    if simBehavior(isubj).rslt(iblock).conf(iparticle,itrial) == 2      % if high confidence
                        temp_signConf(temp_binIndex) = temp_signConf(temp_binIndex) + 1;% add to appopriate signed llr bin
                    end
                end
                
            end %loop over trials
            
        end %loop over blocks
        temp_signConf./temp_binTots;
        % update holding structures
        if condition == 1
            modelStruct(isubj).sel.past.rev         = sum(tempRevRate)./size(tempRevRate,1); %1
            modelStruct(isubj).sel.past.revconf     = sum(tempRevConf)./size(tempRevConf,1); %1
            modelStruct(isubj).sel.past.rep         = temp_binReps./temp_binTots; %2
            modelStruct(isubj).sel.past.signconf    = temp_signConf./temp_binTots; %2 high confidence rate on signed llr
        elseif condition == 2
            modelStruct(isubj).sel.future.rev       = sum(tempRevRate)./size(tempRevRate,1); %1
            modelStruct(isubj).sel.future.revconf   = sum(tempRevConf)./size(tempRevConf,1); %1
            modelStruct(isubj).sel.future.rep       = temp_binReps./temp_binTots; %2
            modelStruct(isubj).sel.future.signconf  = temp_signConf./temp_binTots; %2 high confidence rate on signed llr
        end
    end
  
end

%% Save to file

save('modelStruct','modelStruct');
%% Local Functions
function binNum = binIndex(signllr)
% This function assigns the bin number depending on the llr value 
% signllr can be a scalar or vector 

    illr = 1;
    for llr = signllr
        % data structure: binStruct: [(-inf,-3), [-3:-2), [-2:-1), [-1:0), [0:1), [1:2), [2:3) [3:inf)]
        if llr < 50
            tempbinNum = 8;
            if llr < 3
                tempbinNum = 7;
                if llr < 2    
                    tempbinNum = 6;
                    if llr < 1
                        tempbinNum = 5;
                        if llr < 0
                            tempbinNum = 4;
                            if llr < -1
                                tempbinNum = 3;
                                if llr < -2
                                    tempbinNum = 2;
                                    if llr < -3
                                        tempbinNum = 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        binNum(illr) = tempbinNum;
        illr = illr+1;
    end
end
