%% Setup
% Static variables (set and forget)
n_subjects  = 30;   % set to number of subjects we are analyzing
excluded    = [21]; % excluded subject numbers

% Dynamic variables (CHANGE VALUES HERE BASED ON ANALYSIS)
blockFilter     = 'direction';  % 'direction' or 'volatility'
befAftTrials    = 4;            % how many trials before and after change point
exclude         = true;         % set to true to exclude shitty results
skipExclComp    = true;         % set to true to not run the exclusion comparison plots
xAxis           = [-4:.1:4];    % establish points of support for function and plotting
if exclude
    excluded = [14 20 21 22 27];
end

% Plot legend text
if strcmpi(blockFilter, 'direction')
    txtCond1 = 'Postdiction <-';
    txtCond2 = 'Prediction  ->';
else
    txtCond1 = 'Low volatility';
    txtCond2 = 'High volatility';
end

% Temporary structures
tempZeros  = zeros(1,befAftTrials*2); % temporary array for convenience in loop

% Subjects' responses and confidence values over blocks categorized by blockFilter (1 or 2)
signedconfsCtr1 = [];       % counts high confidence for bin in a given block on condition 1
signedconfsCtr2 = [];       % counts high confidence for bin in a given block on condition 2
signedconfsTotCtr1 = [];    % counts trials falling in bin in a given block on condition 1
signedconfsTotCtr2 = [];    % counts trials falling in bin in a given block on condition 2

% Sequence LLRs for all trials organized according to signedconfsCtr layout
doubleSigStruct = struct;
doubleSigStruct.cond1 = struct;
doubleSigStruct.cond2 = struct;

%% Data Filtering Loop
% loop through subjects and extract appropriate data into local memory

% Pointers
subjIndexTracker = 1; % starting index on resps or confs on where a subject begins
for isubj = 1:n_subjects
    
    if ismember(isubj, excluded)
       continue; 
    end
    
    % signed llr and confidence values on positive signed llr
    doubleSigStruct.cond1(isubj).sllr_pos = [];
    doubleSigStruct.cond2(isubj).sllr_pos = [];
    doubleSigStruct.cond1(isubj).conf_pos = [];
    doubleSigStruct.cond2(isubj).conf_pos = [];
    
    % signed llr and confidence values on negative signed llr
    doubleSigStruct.cond1(isubj).sllr_neg = [];
    doubleSigStruct.cond2(isubj).sllr_neg = [];
    doubleSigStruct.cond1(isubj).conf_neg = [];
    doubleSigStruct.cond2(isubj).conf_neg = [];
    
    % repetition tracking on signed llr
    doubleSigStruct.cond1(isubj).repe_pos = [];
    doubleSigStruct.cond2(isubj).repe_pos = [];
    doubleSigStruct.cond1(isubj).repe_neg = [];
    doubleSigStruct.cond2(isubj).repe_neg = [];
    
    %load subject information
    subject  = isubj;
    filename = dir(sprintf('Data/TURFU_S%02d_*.mat', subject));
    filename = filename.name;
    load(sprintf('Data/%s',filename));
    
    for iblock = 3:10 % loop through the real blocks

        tempSignedConfTot   = tempZeros; % temp array to hold counts of all confidence based on signed seqllr
        tempSignedConf      = tempZeros; % temp array to hold counts of HIGH confidence based on signed seqllr
        tempsignedllrs_pos  = [];
        tempconfidence_pos  = [];
        tempsignedllrs_neg  = [];
        tempconfidence_neg  = [];
        temprepeats_pos     = [];
        temprepeats_neg     = [];
        
        % loop through the trials 
        for i = 2:72    
            
            % REPETITIVE CHOICE signed by EVIDENCE 
            % Note: 1-posLLR, 2-negLLR
            %
            % data structure: binStruct: [(-inf,-3), [-3:-2), [-2:-1), [-1:0), [0:1), [1:2), [2:3) [3:inf)]
            % counters      : [repeat] = 1x8 int array
            %                 [total]  = 1x8 int array
            % Algo:
            % for all trials
            %   sign seqllr by previous choice
            %   add to appropriate bin in [total]
            %   if choice(t) == choice(t-1)
            %       add to appropriate bin in [repeat]
            prevSubjChoice   = (expe.rslt(iblock).resp(i)-1.5)*-2; % converts/maps previous choices (1,2) to (+1,-1)
            signedllr = prevSubjChoice*expe.blck(iblock).seqllr(i);
            tempRepIndex = binIndex(signedllr); % find the appropriate bin
            if signedllr < 0
                tempsignedllrs_neg = horzcat(tempsignedllrs_neg, signedllr);
                tempconfidence_neg = horzcat(tempconfidence_neg, expe.rslt(iblock).conf(i+1));
                if expe.rslt(iblock).resp(i+1) == expe.rslt(iblock).resp(i)
                    temprepeats_neg = horzcat(temprepeats_neg, 1);
                else
                    temprepeats_neg = horzcat(temprepeats_neg, 0);
                end
                
            else
                tempsignedllrs_pos = horzcat(tempsignedllrs_pos, signedllr);
                tempconfidence_pos = horzcat(tempconfidence_pos, expe.rslt(iblock).conf(i+1));
                if expe.rslt(iblock).resp(i+1) == expe.rslt(iblock).resp(i)
                    temprepeats_pos = horzcat(temprepeats_pos, 1);
                else
                    temprepeats_pos = horzcat(temprepeats_pos, 0);
                end
            end
            
            % confidence analysis on signed seqllr
            tempSignedConfTot(tempRepIndex) = tempSignedConfTot(tempRepIndex) + 1;  % count in general
            if expe.rslt(iblock).conf(i+1) == 2                                       % if confidence is high
                tempSignedConf(tempRepIndex) = tempSignedConf(tempRepIndex) + 1;    % count if high confidence
            end
        end
        
        % fill in actual holding structures w/ temp arrays by filtering condition
        if strcmpi(blockFilter,'direction')
            if expe.blck(iblock).taskid == 1
                signedconfsCtr1     = vertcat(signedconfsCtr1, tempSignedConf);
                signedconfsTotCtr1  = vertcat(signedconfsTotCtr1, tempSignedConfTot);
                doubleSigStruct.cond1(isubj).sllr_pos = horzcat(doubleSigStruct.cond1(isubj).sllr_pos, tempsignedllrs_pos);
                doubleSigStruct.cond1(isubj).conf_pos = horzcat(doubleSigStruct.cond1(isubj).conf_pos, tempconfidence_pos);
                doubleSigStruct.cond1(isubj).sllr_neg = horzcat(doubleSigStruct.cond1(isubj).sllr_neg, tempsignedllrs_neg);
                doubleSigStruct.cond1(isubj).conf_neg = horzcat(doubleSigStruct.cond1(isubj).conf_neg, tempconfidence_neg);
                doubleSigStruct.cond1(isubj).repe_pos = horzcat(doubleSigStruct.cond1(isubj).repe_pos, temprepeats_pos);
                doubleSigStruct.cond1(isubj).repe_neg = horzcat(doubleSigStruct.cond1(isubj).repe_neg, temprepeats_neg);
            elseif expe.blck(iblock).taskid == 2
                signedconfsCtr2     = vertcat(signedconfsCtr2, tempSignedConf);
                signedconfsTotCtr2  = vertcat(signedconfsTotCtr2, tempSignedConfTot);
                doubleSigStruct.cond2(isubj).sllr_pos = horzcat(doubleSigStruct.cond2(isubj).sllr_pos, tempsignedllrs_pos);
                doubleSigStruct.cond2(isubj).conf_pos = horzcat(doubleSigStruct.cond2(isubj).conf_pos, tempconfidence_pos);
                doubleSigStruct.cond2(isubj).sllr_neg = horzcat(doubleSigStruct.cond2(isubj).sllr_neg, tempsignedllrs_neg);
                doubleSigStruct.cond2(isubj).conf_neg = horzcat(doubleSigStruct.cond2(isubj).conf_neg, tempconfidence_neg);
                doubleSigStruct.cond2(isubj).repe_pos = horzcat(doubleSigStruct.cond2(isubj).repe_pos, temprepeats_pos);
                doubleSigStruct.cond2(isubj).repe_neg = horzcat(doubleSigStruct.cond2(isubj).repe_neg, temprepeats_neg);
            else
                error(['unexpected value of taskid on block ' num2str(iblock)]);
            end
        elseif strcmpi(blockFilter, 'volatility')
            if expe.blck(iblock).condtn == 1
                signedconfsCtr1     = vertcat(signedconfsCtr1, tempSignedConf);
                signedconfsTotCtr1  = vertcat(signedconfsTotCtr1, tempSignedConfTot);
            elseif expe.blck(iblock).condtn == 2
                signedconfsCtr2     = vertcat(signedconfsCtr2, tempSignedConf);
                signedconfsTotCtr2  = vertcat(signedconfsTotCtr2, tempSignedConfTot);
            else
                error(['unexpected value of condtn on block ' num2str(iblock)]);
            end
        else
            error('blockFilter must be set to an appropriate string value.');
        end
    end
    
    % track where in the pile starts any specific subject for single-subject analyses
    if isubj ~= n_subjects
        subjIndexTracker = vertcat(subjIndexTracker, size(signedconfsCtr1,1)+1);
    end
end

%% Percentage of repeated choices over seqLLR signed by choice
% calculate rate of repeated choices from counters in main loop

% means
signedconfrate1 = zeros(1,8);
signedconfrate2 = zeros(1,8);

for i = 1:8
    % mean calculations
    signedconfrate1(i) = sum(signedconfsCtr1(:,i))/sum(signedconfsTotCtr1(:,i));
    signedconfrate2(i) = sum(signedconfsCtr2(:,i))/sum(signedconfsTotCtr2(:,i));
end

% The following code is not optimized for text space. It is separated into redundant
% loops for clarity.

% data organized to generate statistics for subject-wide variability on CONFIDENCE
sssignconfrate1 = zeros(1,8); % single subject rate of high confidence on condition 1
sssignconfrate2 = zeros(1,8); % single subject rate of high confidence on condition 2
subjCtr = 1;
for i = 1:numel(signedconfsCtr1(:,1))
        sssignconfrate1(subjCtr,:) = sssignconfrate1(subjCtr,:) + signedconfsCtr1(i,:);
        sssignconfrate2(subjCtr,:) = sssignconfrate2(subjCtr,:) + signedconfsCtr2(i,:);
    if mod(i,4) == 0
        sssignconfrate1(subjCtr,:) = sssignconfrate1(subjCtr,:)./sum(signedconfsTotCtr1(i-3:i,:));
        sssignconfrate2(subjCtr,:) = sssignconfrate2(subjCtr,:)./sum(signedconfsTotCtr2(i-3:i,:));
        if i ~= numel(signedconfsCtr1(:,1))
            sssignconfrate1 = vertcat(sssignconfrate1,zeros(1,8));
            sssignconfrate2 = vertcat(sssignconfrate2,zeros(1,8));
        end
        subjCtr = subjCtr + 1;
    end
end

%% statistics on confidence in the two conditions
for i = 1:size(sssignconfrate1,2)
    [pConfSigned(i) hConfSigned(i)  statsConfSigned(i)] = signrank(sssignconfrate1(:,i), sssignconfrate2(:,i));
end
statsConfSignedAsterisks = {};
statCtr = 1;
for stat = pConfSigned
    if stat < .001 
        statsConfSignedAsterisks = [statsConfSignedAsterisks '***'];
    elseif stat < .01
        statsConfSignedAsterisks = [statsConfSignedAsterisks '**'];
    elseif stat < .05
        statsConfSignedAsterisks = [statsConfSignedAsterisks '*'];
    else
        statsConfSignedAsterisks = [statsConfSignedAsterisks ' '];
    end
    statCtr = statCtr + 1;
end


%% Fit sigmoid on p(stay) with signed llr as predictor 
% this is to be used to weigh the two confidence sigmoids
%
% Logistic curve parameters for all subjects, 2 conditions
paramsrepecond1 = zeros(n_subjects,2); 
paramsrepecond2 = zeros(n_subjects,2);

% discrete support for the curve (for plotting later)
repeatcurve_cond1 = zeros(n_subjects,numel(xAxis));
repeatcurve_cond2 = zeros(n_subjects,numel(xAxis));

for isubj = 1:n_subjects
    
    if ismember(isubj, excluded)
       continue; 
    end
    x_cond1 = [];
    n_cond1 = [];
    y_cond1 = [];
    x_cond2 = [];
    n_cond2 = [];
    y_cond2 = [];
    
    % concatenate the separated positive and negative llr data on repeats together
    x_cond1 = [doubleSigStruct.cond1(isubj).sllr_pos doubleSigStruct.cond1(isubj).sllr_neg];
    n_cond1 = [ones(1,numel(doubleSigStruct.cond1(isubj).repe_pos)) ones(1,numel(doubleSigStruct.cond1(isubj).repe_neg))];
    y_cond1 = [doubleSigStruct.cond1(isubj).repe_pos doubleSigStruct.cond1(isubj).repe_neg];
    
    x_cond2 = [doubleSigStruct.cond2(isubj).sllr_pos doubleSigStruct.cond2(isubj).sllr_neg];
    n_cond2 = [ones(1,numel(doubleSigStruct.cond2(isubj).repe_pos)) ones(1,numel(doubleSigStruct.cond2(isubj).repe_neg))];
    y_cond2 = [doubleSigStruct.cond2(isubj).repe_pos doubleSigStruct.cond2(isubj).repe_neg];
    
    paramsrepecond1(isubj,:) = glmfit(x_cond1, [y_cond1.' n_cond1.'],'binomial','link','logit');
    paramsrepecond2(isubj,:) = glmfit(x_cond2, [y_cond2.' n_cond2.'],'binomial','link','logit');
    
    z_cond1 = paramsrepecond1(isubj,1)+(paramsrepecond1(isubj,2)*xAxis);
    z_cond1 = 1 ./ (1 + exp(-z_cond1));
    repeatcurve_cond1(isubj,:) = z_cond1;
    
    z_cond2 = paramsrepecond2(isubj,1)+(paramsrepecond2(isubj,2)*xAxis);
    z_cond2 = 1 ./ (1 + exp(-z_cond2));
    repeatcurve_cond2(isubj,:) = z_cond2;
    
end

%% Fit two sigmoids on confidence choices with signed llrs as predictors on either side of 0
%   One sigmoid is fitted to confidence values on evidence to stay
%   Other sigmoid is fitted to confidence values on evidence to switch

% hold parameters for the two conditions on the two sides of 0
paramsconfcond1_pos = zeros(n_subjects,2);
paramsconfcond1_neg = zeros(n_subjects,2);
paramsconfcond2_pos = zeros(n_subjects,2);
paramsconfcond2_neg = zeros(n_subjects,2);

% the double sigmoid curves on confidence, separated
confcurve_cond1_pos = zeros(n_subjects,numel(xAxis));
confcurve_cond2_pos = zeros(n_subjects,numel(xAxis));
confcurve_cond1_neg = zeros(n_subjects,numel(xAxis));
confcurve_cond2_neg = zeros(n_subjects,numel(xAxis));

for isubj = 1:n_subjects
    
    if ismember(isubj, excluded)
       continue; 
    end
    % fit over evidence to stay, condition 1
    x = doubleSigStruct.cond1(isubj).sllr_pos;
    n = ones(1,numel(doubleSigStruct.cond1(isubj).conf_pos));
    y = doubleSigStruct.cond1(isubj).conf_pos - 1;
    paramsconfcond1_pos(isubj,:) = glmfit(x, [y.' n.'],'binomial','link','logit');
    z = paramsconfcond1_pos(isubj,1)+(paramsconfcond1_pos(isubj,2)*xAxis);
    z = 1 ./ (1 + exp(-z));
    confcurve_cond1_pos(isubj,:) = z;
    
    % fit over evidence to switch, condition 1
    x = doubleSigStruct.cond1(isubj).sllr_neg;
    n = ones(1,numel(doubleSigStruct.cond1(isubj).conf_neg));
    y = doubleSigStruct.cond1(isubj).conf_neg - 1;
    paramsconfcond1_neg(isubj,:) = glmfit(x, [y.' n.'],'binomial','link','logit');
    z = paramsconfcond1_neg(isubj,1)+(paramsconfcond1_neg(isubj,2)*xAxis);
    z = 1 ./ (1 + exp(-z));
    confcurve_cond1_neg(isubj,:) = z;
    
    % fit over evidence to stay, condition 2
    x = doubleSigStruct.cond2(isubj).sllr_pos;
    n = ones(1,numel(doubleSigStruct.cond2(isubj).conf_pos));
    y = doubleSigStruct.cond2(isubj).conf_pos - 1;
    paramsconfcond2_pos(isubj,:) = glmfit(x, [y.' n.'],'binomial','link','logit');
    z = paramsconfcond2_pos(isubj,1)+(paramsconfcond2_pos(isubj,2)*xAxis);
    z = 1 ./ (1 + exp(-z));
    confcurve_cond2_pos(isubj,:) = z;
    
    % fit over evidence to switch, condition 2
    x = doubleSigStruct.cond2(isubj).sllr_neg;
    n = ones(1,numel(doubleSigStruct.cond2(isubj).conf_neg));
    y = doubleSigStruct.cond2(isubj).conf_neg - 1;
    paramsconfcond2_neg(isubj,:) = glmfit(x, [y.' n.'],'binomial','link','logit');
    z = paramsconfcond2_neg(isubj,1)+(paramsconfcond2_neg(isubj,2)*xAxis);
    z = 1 ./ (1 + exp(-z));
    confcurve_cond2_neg(isubj,:) = z;
end

%% Find cutoff values for low/high confidence (SKIP if not trying to do this)
if false
% Cutoff values in the two conditions
confCutoff = [];
confCutoff(1,:,:) = zeros(n_subjects,2); % cutoffs for condition 1
confCutoff(2,:,:) = zeros(n_subjects,2); % cutoffs for condition 2
%BADs settings
addpath(genpath('bads-master')); % load the BADS files into working directory
options = [];
options.TolMesh = 0.1;
for isubj = 1:n_subjects
    if ismember(isubj, excluded)
       continue; 
    end
    for icond = 1:2
        % make parameters easier to access for function
        if icond == 1
            param_stay = paramsrepecond1(isubj,:);
            param_neg  = paramsconfcond1_neg(isubj,:);
            param_pos  = paramsconfcond1_pos(isubj,:);
        else
            param_stay = paramsrepecond2(isubj,:);
            param_neg  = paramsconfcond2_neg(isubj,:);
            param_pos  = paramsconfcond2_pos(isubj,:);
        end
        % store cutoff value for p(high confidence) = .5 on negative llr
        
        %BADS settings
        startPt = [rand()-1]; % vector corresponding to values of parameters to be estimated
        lBound  = [-7]; % HARD lower bound of parameter values
        uBound  = [0]; % HARD upper bound of parameter values
        pLBound = [-7]; % Plausible lower bound of parameter values 
        pUBound = [0]; % Plausible upper bound of parameter values

        [optiParams, fVal] = bads(@cutoffValue, startPt, lBound, uBound, pLBound, pUBound, [], options);
        confCutoff(icond,isubj,1) = optiParams;

        % store cutoff value for p(high confidence) = .5 on positive llr
        
        %BADS settings
        startPt = [rand()]; % vector corresponding to values of parameters to be estimated
        lBound  = [0]; % HARD lower bound of parameter values
        uBound  = [7]; % HARD upper bound of parameter values
        pLBound = [0]; % Plausible lower bound of parameter values (optional)
        pUBound = [7]; % Plausible upper bound of parameter values (optional)

        [optiParams, fVal] = bads(@cutoffValue, startPt, lBound, uBound, pLBound, pUBound, [], options);
        confCutoff(icond,isubj,2) = optiParams;
    end 
end

%sanity check
disp(['Subject ' ' Condition 1 (lo/hi) ' ' Condition 2 (lo/hi) ']);
for isubj = 1:n_subjects
    if ismember(isubj, excluded)
       continue; 
    end
    disp([num2str(isubj) ' ' num2str(confCutoff(1,isubj,1)) ' ' num2str(confCutoff(1,isubj,2)) ' '...
          num2str(confCutoff(2,isubj,1)) ' ' num2str(confCutoff(2,isubj,2))]);
end
%save('confCutoff','confCutoff');
end

%% condition 1 (postdiction)
disp('postdiction smoothed double sigmoid confidence thing');
collect_curves_cond1 = [];
avgCurve_cond1 = [];
semCurve_cond1 = [];
for isubj = 1:n_subjects
    if ismember(isubj, excluded)
       continue; 
    end
    z = dblSmoothSigmoid(confcurve_cond1_neg(isubj,:),confcurve_cond1_pos(isubj,:),repeatcurve_cond1(isubj,:));
    collect_curves_cond1 = vertcat(collect_curves_cond1, z);
    % to check single subject confidence curves
    %{ 
    figure(isubj);
    hold on;
    ylim([0 1]);
    yline(.5);
    plot(xAxis, z);
    hold off;
    pause(.5);
    %}
end
avgCurve_cond1 = mean(collect_curves_cond1);
semCurve_cond1 = std(collect_curves_cond1)./sqrt(n_subjects-numel(excluded));

%% condition 2 (prediction)
disp('prediction smoothed double sigmoid confidence thing');
collect_curves_cond2 = [];
avgCurve_cond2 = [];
semCurve_cond2 = [];
for isubj = 1:n_subjects
    if ismember(isubj, excluded)
       continue; 
    end
    z = dblSmoothSigmoid(confcurve_cond2_neg(isubj,:),confcurve_cond2_pos(isubj,:),repeatcurve_cond2(isubj,:));
    collect_curves_cond2 = vertcat(collect_curves_cond2, z);
    % to check single subject confidence curves
    %{ 
    figure(isubj);
    hold on;
    ylim([0 1]);
    yline(.5);
    plot(xAxis, z);
    hold off;
    pause(.5);
    %}
end
avgCurve_cond2 = mean(collect_curves_cond2);
semCurve_cond2 = std(collect_curves_cond2)./sqrt(n_subjects-numel(excluded));

%% Parameter comparison

% remove 0's from data structs
slopes_neg1 = paramsconfcond1_neg(:,2);
slopes_neg2 = paramsconfcond2_neg(:,2);
slopes_pos1 = paramsconfcond1_pos(:,2);
slopes_pos2 = paramsconfcond2_pos(:,2);

slopes_neg1(slopes_neg1==0)=[];
slopes_neg2(slopes_neg2==0)=[];
slopes_pos1(slopes_pos1==0)=[];
slopes_pos2(slopes_pos2==0)=[];

disp('comparing slopes on negative evidence side');
[h p stats] = signrank(slopes_neg1, slopes_neg2)
figure(1);
hold on;
negMeanPlot     = bar([1:2], [mean(slopes_neg1) mean(slopes_neg2)]);
negErrorPlot    = errorbar([1:2], [mean(slopes_neg1) mean(slopes_neg2)], [std(slopes_neg1) std(slopes_neg2)]./sqrt(25),...
                           'LineStyle','none','CapSize',0);                    
xticks([1:2]);
xticklabels({'postdiction' 'prediction'});
ylim([min(mean(slopes_neg1),mean(slopes_neg2))-.1 0]);   
hold off;

disp('comparing slopes on positive evidence side');
[h p stats] = signrank(slopes_pos1, slopes_pos2)
figure(2);
hold on;
posMeanPlot     = bar([1:2], [mean(slopes_pos1) mean(slopes_pos2)]);
posErrorPlot    = errorbar([1:2], [mean(slopes_pos1) mean(slopes_pos2)], [std(slopes_pos1) std(slopes_pos2)]./sqrt(25),...
                           'LineStyle','none','CapSize',0);
xticks([1:2]);
xticklabels({'postdiction' 'prediction'});
ylim([0 max(mean(slopes_pos1),mean(slopes_pos2))+.1]);   
hold off;

%% Plotting
figure(30);
hold on;
% condition 1 stdev over single subject means
p30_1s = errorbar([1:8], mean(sssignconfrate1),std(sssignconfrate1)/sqrt(length(sssignconfrate1)),'LineStyle','none',...
                'Color',[1 .75 .5],'LineWidth',1);    
% condition 1 mean over single subject means
p30_1a = scatter([1:8], mean(sssignconfrate1),'filled','MarkerFaceColor',[1 .75 .5],'MarkerEdgeColor',...
                 [1 .75 .5]);  
% condition 2 stdev over single subject means
p30_2s = errorbar([1:8], mean(sssignconfrate2),std(sssignconfrate2)/sqrt(length(sssignconfrate2)),'LineStyle','none',...
                'Color',[.75 .75 1],'LineWidth',1);   
% condition 2 mean over single subject means
p30_2a = scatter([1:8], mean(sssignconfrate2),'filled','d','MarkerFaceColor',[.75 .75 1],'MarkerEdgeColor',...
                 [.75 .75 1]); 

asterisksY = mean(sssignconfrate1);     % statistical significance over interval
for i = 1:numel(statsConfSignedAsterisks)
    text(i, asterisksY(i)+.08, cell2mat(statsConfSignedAsterisks(i)),'FontSize',14);
end
%title({sprintf('Effect of %s on confidence vs evidence signed by the previous choice', blockFilter) '25 subjects; error bars SEM'});
xticks([1:8]);
xticklabels({'(-\infty,-3)', '[-3:-2)', '[-2:-1)', '[-1:0)', '[0:1)', '[1:2)', '[2:3)', '[3:\infty)'}); 
xlabel('Strength of evidence (LLR) at trial t for(+) or against(-) subject prior belief at trial t-1');
ylabel('Proportion of trials with high subject confidence');
xlim([1 8]);
ylim([.2 1]);

% plot smoothed double sigmoid curves for both conditions
model_1 = plot(xAxis+4.5, avgCurve_cond1,'LineStyle','--','Color',[1 .75 .5]);
model_2 = plot(xAxis+4.5, avgCurve_cond2,'LineStyle','--','Color',[.75 .75 1]);
err1 = shadedErrorBar(xAxis+4.5, avgCurve_cond1, semCurve_cond1,'lineprops',{'--','Color',[1 .75 .5]},'transparent',1);
err2 = shadedErrorBar(xAxis+4.5, avgCurve_cond2, semCurve_cond2,'lineprops',{'--','Color',[.75 .75 1]},'transparent',1);

legend([p30_1a p30_2a model_1 model_2], {'data postdictive','data predictive','model postdictive','model predictive'},'Location', 'southeast')
set(gca, 'FontName', 'Helvetica');
hold off;

%% Local Functions
function y = dblSmoothSigmoid(sgmd_neg, sgmd_pos, repeat)
    % This function weighs the two sigmoid curves on the confidence against and for prior
    % choice by the proportion of repetition when that confidence was made
    y = sgmd_neg.*(1-repeat) + sgmd_pos.*repeat;
end
    
function objFn = cutoffValue(llr)
    param_stay = evalin('base', 'param_stay');
    param_neg  = evalin('base', 'param_neg');
    param_pos  = evalin('base', 'param_pos');
    
    pstay = param_stay(1)+param_stay(2).*llr;
    pstay = 1./(1+exp(-pstay));
    
    lneg = param_neg(1)+param_neg(2).*llr;
    lneg = 1./(1+exp(-lneg));
    
    lpos = param_pos(1)+param_pos(2).*llr;
    lpos = 1./(1+exp(-lpos));
    
    cutoff = (lneg.*(1-pstay))+(lpos.*pstay);
    
    objFn = (.5-cutoff)^2;
end

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
