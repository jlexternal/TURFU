% Convert model posterior beliefs to confidence values
%
%   Name:   model_belief2conf.m
%   Type:   Script
% Author:   Jun Seok Lee
%   Date:   May 2019
% 

% In order for this code to function, the following structures are needed:
% 1. simBehavior structures from 'model'_sim.m should have already been run and stored
%    into savedStructure folder 
%    ('model' a placeholder for the name of model)
% 2. confCutoff matrix from conf_doubleSigmoid.m 

% This code will add to the simBeehavior structure for each model another field
% corresponding to the model's confidence at each trial


%% Setup

% Static variables (set and forget)
n_subjects  = 30;   % set to number of subjects we are analyzing
excluded    = [21]; % excluded subject numbers

% **Dynamic variables (CHANGE VALUES HERE BASED ON ANALYSIS)**
blockFilter     = 'direction'; % 'direction' or 'volatility'
condition       = 2;           % 1 for past/low; 2 for future/high (direction/volatility)
befAftTrials    = 4;           % how many trials before and after change point
exclude         = true;        % set to true to exclude shitty results
if exclude
    excluded = [14 20 21 22 27];
end

%load saved confCutoff matrix 
filename = dir('SavedStructures/confCutoff.mat');
filename = filename.name;
load(sprintf('SavedStructures/%s',filename));

%names of the simBehavior structures for the three models
if condition == 1
    modelFileStrings = ["glazeOriginal_sim_post_struct_*.mat","infNoiseGlaze_sim_post_struct_*.mat","selNoiseGlaze_sim_post_struct_*.mat"];
elseif condition == 2
    modelFileStrings = ["glazeOriginal_sim_pred_struct_*.mat","infNoiseGlaze_sim_pred_struct_*.mat","selNoiseGlaze_sim_pred_struct_*.mat"];
end

for modelStr = modelFileStrings
    %load the appropriate simBehavior from each model
    filename = dir(strcat('SavedStructures/',modelStr));
    filename = filename.name;
    load(sprintf('SavedStructures/%s',filename));
    
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
        %load cutoff values for subject 
        cutoff_neg = confCutoff(condition, isubj,1);
        cutoff_pos = confCutoff(condition, isubj,2);

        for iblock = 3:10 % go through the blocks corresponding to the condition of interest
            if expe.blck(iblock).taskid ~= condition 
                continue; 
            end
            % go through the beliefs generated by the model and convert to confidence
            beliefs = simBehavior(isubj).rslt(iblock).belief;
            prevChoices = simBehavior(isubj).rslt(iblock).resp;
            
            simBehavior(isubj).rslt(iblock).conf = conv2conf(beliefs, prevChoices, cutoff_neg, cutoff_pos);
        end
    end
    savename = dir(strcat('SavedStructures/',modelStr));
    savename = savename.name;
    save(strcat('SavedStructures/',savename),'simBehavior');
end

%% Local Functions

function conf = conv2conf(beliefs, prevChoices, negCut, posCut)
% beliefs is an array of LLRs on a given block
% prevChoices is an array of responses (1 or 2) on a given block
% negCut is the negative llr cutoff for high confidence found from the weighted double sigmoid
% posCut is the positive llr cutoff for high confidence found from the weighted double sigmoid
    n_particles = numel(beliefs(:,1));
    for i = 1:numel(beliefs(1,:))   
        if i == 1
            % base 1st confidence on positive evidence cutoff
            A = zeros(n_particles,1);
            B = zeros(n_particles,1);
            A(abs(beliefs(:,1)) > posCut) = 2;
            B(abs(beliefs(:,1)) <= posCut) = 1;
            conf(:,1) = A+B;
        else
            % everything else is based on signed (by previous choice) evidence
            Aneg = zeros(n_particles,1);
            Apos = zeros(n_particles,1);
            
            prevsign = (prevChoices(:,i)-1.5).*-2; 
            belief = beliefs(:,i).*prevsign; % belief signed by previous choice
            
            negCutArr = belief<0;
            negbelief = belief.*negCutArr;
            
            posCutArr = belief>0;
            posbelief = belief.*posCutArr;
            
            Aneg(negbelief<=negCut) = 2;
            Aneg(negbelief < 0 & negbelief>negCut)  = 1;
            if posCut ~= 0 
                Apos(posbelief>=posCut) = 2;
                Apos(posbelief > 0 & posbelief<posCut)  = 1;
            else
                Apos(posbelief>posCut) = 2;
            end
            
            conf(:,i) = Aneg+Apos;
        end
    end
end