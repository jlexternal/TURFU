% Code to measure model fit on data based on the 3 competing models. 
% Initial data stored in file 'tempstuff.mat'

clearvars -except infNoiseMeanRespRate infNoiseInclRespRate glazeOrigMeanRespRate glazeOrigInclRespRate selNoiseMeanRespRate selNoiseInclRespRate;
% Static variables (set and forget)
n_subjects  = 30;   % set to number of subjects we are analyzing
excluded    = [21]; % excluded subject numbers

% Dynamic variables (CHANGE VALUES HERE BASED ON ANALYSIS)
blockFilter     = 'direction';  % 'direction' or 'volatility'
exclude         = true;         % set to true to exclude shitty results
specSubjects    = false;        % set to true if analyzing only certain subjects
condition      = [1];          % set to 1 if for 1st condition, 2 for 2nd, 1:2 for both
nparticles = 5000;

if condition == 1
    modelFileStrings = ["glazeOriginal_sim_post_struct_*.mat","infNoiseGlaze_sim_post_struct_*.mat","selNoiseGlaze_sim_post_struct_*.mat"];
elseif condition == 2
    modelFileStrings = ["glazeOriginal_sim_pred_struct_*.mat","infNoiseGlaze_sim_pred_struct_*.mat","selNoiseGlaze_sim_pred_struct_*.mat"];
end

glazeorigchoice1    = zeros(30,8);
glazeorigchoicetot	= zeros(30,8);
infnoisechoice1   	= zeros(30,8);
infnoisechoicetot   = zeros(30,8);
selnoisechoice1     = zeros(30,8);
selnoisechoicetot   = zeros(30,8);

%% Organize simulated data

for modelStr = modelFileStrings
    %load the appropriate simBehavior from each model
    filename = dir(strcat('SavedStructures/',modelStr));
    filename = filename.name;
    load(sprintf('SavedStructures/%s',filename));
    
    for cond = condition
        blocks = find([expe.blck.condtn] == cond);
        for block = blocks
            
            
            for itrial = 1:72
                % update number of total particles for the seqllr bin on the trial
                if strcmpi(modelStr,modelFileStrings(1))
                    tempResp = simBehavior(isubj).rslt(block).resp(itrial+1);
                    tempChoicetot = glazeorigchoicetot(isubj,binIndex(expe.blck(block).seqllr(itrial)));
                    glazeorigchoicetot(isubj,binIndex(expe.blck(block).seqllr(itrial))) = tempChoicetot + 1;
                    if tempResp == 1
                        tempChoice1 = glazeorigchoice1(isubj,binIndex(expe.blck(block).seqllr(itrial)));
                        glazeorigchoice1(isubj,binIndex(expe.blck(block).seqllr(itrial))) = tempChoice1 + 1;
                    end
                end
                if strcmpi(modelStr,modelFileStrings(2))
                    tempChoicetot = infnoisechoicetot(isubj,binIndex(expe.blck(block).seqllr(itrial)));
                    infnoisechoicetot(isubj,binIndex(expe.blck(block).seqllr(itrial))) = tempChoicetot + nparticles;

                    tempChoice1 = infnoisechoice1(isubj,binIndex(expe.blck(block).seqllr(itrial))); % get previous value in bin
                    tempTrialChoice1 = simBehavior(isubj).rslt(block).resp(:,itrial+1);     % extract array from structure for easiness
                    infnoisechoice1(isubj,binIndex(expe.blck(block).seqllr(itrial))) = tempChoice1 + numel(tempTrialChoice1(tempTrialChoice1==1)); % update bin
                end
                if strcmpi(modelStr,modelFileStrings(3))
                    tempChoicetot = selnoisechoicetot(isubj,binIndex(expe.blck(block).seqllr(itrial)));
                    selnoisechoicetot(isubj,binIndex(expe.blck(block).seqllr(itrial))) = tempChoicetot + nparticles;
                    tempChoice1 = selnoisechoice1(isubj,binIndex(expe.blck(block).seqllr(itrial))); % get previous value in bin
                    tempTrialChoice1 = simBehavior(isubj).rslt(block).resp(:,itrial+1);     % extract array from structure for easiness
                    selnoisechoice1(isubj,binIndex(expe.blck(block).seqllr(itrial))) = tempChoice1 + numel(tempTrialChoice1(tempTrialChoice1==1)); % update bin
                end
                
            end
        end
    end
end


%% Organize subject data
if exclude
    excluded = [14 20 21 22 27];
end
if specSubjects
    subjects = [1];          % specify subjects here
else
    subjects = [1:n_subjects];
end

subjdatachoice1     = zeros(30,8);
subjdatachoicetot   = zeros(30,8);

for isubj = 1:n_subjects
    
    if ismember(isubj, excluded)
       continue; 
    end
    
    %load subject information
    subject  = isubj;
    filename = dir(sprintf('Data/TURFU_S%02d_*.mat', subject));
    filename = filename.name;
    load(sprintf('Data/%s',filename));
    
    for cond = condition
        blocks = find([expe.blck.condtn] == cond);
        for block = blocks
            
            % get p(choice1) curves on seqllr
            for itrial = 1:72
                tempResp = expe.rslt(block).resp(itrial+1);
                tempChoicetot = subjdatachoicetot(isubj,binIndex(expe.blck(block).seqllr(itrial)));
                subjdatachoicetot(isubj,binIndex(expe.blck(block).seqllr(itrial))) = tempChoicetot + 1;
                if tempResp == 1
                    tempChoice1 = subjdatachoice1(isubj,binIndex(expe.blck(block).seqllr(itrial)));
                    subjdatachoice1(isubj,binIndex(expe.blck(block).seqllr(itrial))) = tempChoice1 + 1;
                end
            end
        end
    end
end
%% Convert to percentages and plot

respRate = subjdatachoice1./subjdatachoicetot;
exclRows = isnan(respRate(:,1));
inclRespRate = respRate(~exclRows,:);
meanRespRate = mean(inclRespRate);

%scatter([1:8], meanRespRate);
hold on;
curve_data = errorbar([1:8], meanRespRate, std(inclRespRate),'k','LineWidth',1.5);
%scatter([1:8], infNoiseMeanRespRate);
curve_glazeOrig = errorbar([1:8], glazeOrigMeanRespRate, std(glazeOrigInclRespRate),'r','LineWidth',1.4);
curve_infNoise  = errorbar([1:8], infNoiseMeanRespRate, std(infNoiseInclRespRate),'b','LineWidth',1.3);
curve_selNoise  = errorbar([1:8], selNoiseMeanRespRate, std(selNoiseInclRespRate),'g','LineWidth',1.2);
legend([curve_data,curve_glazeOrig,curve_infNoise,curve_selNoise],...
       {'Data','Deterministic','+infNoise','+infNoise+selNoise'},'Location','southeast');
ylim([0 1]);
xticks([1:8]);
xticklabels({'(-\infty,-3)', '[-3:-2)', '[-2:-1)', '[-1:0)', '[0:1)', '[1:2)', '[2:3)', '[3:\infty)'}); 
xlabel('Sequence LLR');
ylabel('Proportion of Choice corresponding to +''ve LLR');
hold off;


%% Local Functions
function binNum = binIndex(signllr)
    % This function assigns the bin number depending on the llr value signed by
    %   the previous subject choice

    % data structure: binStruct: [(-inf,-3), [-3:-2), [-2:-1), [-1:0), [0:1), [1:2), [2:3) [3:inf)]
    if signllr < 50
        binNum = 8;
        if signllr < 3
            binNum = 7;
            if signllr < 2    
                binNum = 6;
                if signllr < 1
                    binNum = 5;
                    if signllr < 0
                        binNum = 4;
                        if signllr < -1
                            binNum = 3;
                            if signllr < -2
                                binNum = 2;
                                if signllr < -3
                                    binNum = 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

    