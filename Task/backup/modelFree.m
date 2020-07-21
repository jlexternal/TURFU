% Model free analysis on TURFU subjects
%
%   Name:   modelFree.m
%   Type:   Script
% Author:   Jun Seok Lee
%   Date:   April 2019
% 
% Description:
%   Puts all responses around the reversal (change) point for all subjects 
%   into one data structure. Specificially, data structures for intra-subject analysis 
%   (single subject) i.e. 8 trials (4 before, 4 after change point)
%   
%   Responses and confidence on signed sequence LLR curves.

% Update: April 30, 2019
%   Can also be used to compare data to model after model simulations performed.
%   Currently, only comparing responses from data to model.
%   Will compare confidence once confidence thresholds are set.
%   
%   Trying to look at interactions, I am also including some really shitty inefficient code.

clf;
close all;
% a whole lotta colors
postRGB = [1 .75 .5];
predRGB = [.75 .75 1];
detRGB = [0, 0.4470, 0.7410];
infRGB = [0.8500, 0.3250, 0.0980];
selRGB = [0.4660, 0.6740, 0.1880];

%% Setup
% Static variables (set and forget)
n_subjects  = 30;   % set to number of subjects we are analyzing
excluded    = [21]; % excluded subject numbers

% Dynamic variables (CHANGE VALUES HERE BASED ON ANALYSIS)
blockFilter     = 'direction'; % 'direction' or 'volatility'
interaction     = false;       % set to true for 2x2 interaction analysis; blockFilter should be set to 'direction' for this
befAftTrials    = 4;           % how many trials before and after change point
exclude         = true;        % set to true to exclude shitty results
skipExclComp    = true;        % set to true to not run the exclusion comparison plots
data2modelComp = true;         % set to true to only check condition 1 data to condition 1 models

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
resps1 = []; 
resps2 = [];
revconfs1 = [];
revconfs2 = [];

repeat1 = [];
repeat2 = [];
reptot1 = [];
reptot2 = [];
repconfsCtr1 = [];
repconfsCtr2 = [];
repconfsTotCtr1 = [];
repconfsTotCtr2 = [];

signedconfsCtr1 = [];       % counts high confidence for bin in a given block on condition 1
signedconfsCtr2 = [];       % counts high confidence for bin in a given block on condition 2
signedconfsTotCtr1 = [];    % counts trials falling in bin in a given block on condition 1
signedconfsTotCtr2 = [];    % counts trials falling in bin in a given block on condition 2

% Directions (correct response) from the experiment over blocks categorized by blockFilter
expResps1 = [];
expResps2 = [];

% Data structures if interactions are considered
if interaction
    %Past-low
    resps1_lo               = [];
    revconfs1_lo            = [];
    expResps1_lo            = [];

    reptot1_lo              = [];
    repeat1_lo              = [];
    repconfsCtr1_lo         = [];
    repconfsTotCtr1_lo      = [];
    signedconfsCtr1_lo      = [];
    signedconfsTotCtr1_lo   = [];
    %Past-high
    resps1_hi               = [];
    revconfs1_hi            = [];
    expResps1_hi            = [];

    reptot1_hi              = [];
    repeat1_hi              = [];
    repconfsCtr1_hi         = [];
    repconfsTotCtr1_hi      = [];
    signedconfsCtr1_hi      = [];
    signedconfsTotCtr1_hi   = [];
    %Future-low
    resps2_lo               = [];
    revconfs2_lo            = [];
    expResps2_lo            = [];

    reptot2_lo              = [];
    repeat2_lo              = [];
    repconfsCtr2_lo         = [];
    repconfsTotCtr2_lo      = [];
    signedconfsCtr2_lo      = [];
    signedconfsTotCtr2_lo   = [];
    %Future-high
    resps2_hi               = [];
    revconfs2_hi            = [];
    expResps2_hi            = [];

    reptot2_hi              = [];
    repeat2_hi              = [];
    repconfsCtr2_hi         = [];
    repconfsTotCtr2_hi      = [];
    signedconfsCtr2_hi      = [];
    signedconfsTotCtr2_hi   = [];
end

%% Data Filtering Loop
% loop through subjects and extract appropriate data into local memory

% Pointers
subjIndexTracker = [1]; % starting index on resps or confs on where a subject begins

for isubj = 1:n_subjects
    
    if ismember(isubj, excluded)
       continue; 
    end
    
    %load subject information
    subject  = isubj;
    filename = dir(sprintf('Data/TURFU_S%02d_*.mat', subject));
    filename = filename.name;
    load(sprintf('Data/%s',filename));
    
    for iblock = 3:10 % loop through the real blocks
        
        tempRevResp = []; % temp array to hold responses around reversal point
        tempRevConf = []; % temp array to hold confidence around reversal point
        expeRevResp = []; % temp array to hold correct responses around reversal point
        
        tempRepeat          = zeros(1,8); % temp array to hold repeat data structure per block
        tempRepTot          = zeros(1,8); % temp array to hold counts of seqllr signed by choice 
        tempRepConfTot      = zeros(1,8); % temp array to hold counts of all confidence based on signed seqllr
        tempRepConf         = zeros(1,8); % temp array to hold counts of HIGH confidence based on signed seqllr
        tempSignedConfTot   = zeros(1,8); % temp array to hold counts of all confidence based on signed seqllr
        tempSignedConf      = zeros(1,8); % temp array to hold counts of HIGH confidence based on signed seqllr
         
        % loop through the trials 
        for i = 3:72    % start from 3 since starting from 2 will compare with the null choice
            
            % REVERSAL CURVE analysis
            % identify trials after switch 
            if expe.blck(iblock).seqdir(i) ~= expe.blck(iblock).seqdir(i-1)
                % get subject response and confidence
                tempRevResp     = vertcat(tempRevResp, expe.rslt(iblock).resp(i-(befAftTrials-1):i+befAftTrials));
                tempRevConf     = vertcat(tempRevConf, expe.rslt(iblock).conf(i-(befAftTrials-1):i+befAftTrials));
                tempZeros(1,:)  = expe.blck(iblock).seqdir(i-1); % put trial direction before reversal
                expeRevResp     = vertcat(expeRevResp, tempZeros); % put this into the experiment direction
            end
            
            % REPETITIVE CHOICE signed by EVIDENCE 
            % Note: 1-posLLR, 2-negLLR
            %
            % data structure: binStruct: [(-inf,-3), [-3:-2), [-2:-1), [-1:0), [0:1), [1:2), [2:3) [3:inf)]
            % counters      : [repeat] = 1x8 int array
            %                 [total]  = 1x8 int array
            % Algo:
            % for all trials
            %   sign seqllr by choice
            %   add to appropriate bin in [total]
            %   if choice(t) == choice(t-1)
            %       add to appropriate bin in [repeat]
            subjChoice   = (expe.rslt(iblock).resp(i-1)-1.5)*-2; % converts/maps choices (1,2) to (+1,-1)
            
            tempRepIndex = binIndex(subjChoice*expe.blck(iblock).seqllr(i-1)); % find the appropriate bin
            tempRepTot(tempRepIndex) = tempRepTot(tempRepIndex) + 1;            % add to that bin in [total]
            
            % confidence analysis on signed seqllr
            tempSignedConfTot(tempRepIndex) = tempSignedConfTot(tempRepIndex) + 1;  % count in general
            if expe.rslt(iblock).conf(i) == 2                                       % if confidence is high
                tempSignedConf(tempRepIndex) = tempSignedConf(tempRepIndex) + 1;    % count if high confidence
            end

            % choice analysis
            if expe.rslt(iblock).resp(i) == expe.rslt(iblock).resp(i-1)         % if choice repeated
                tempRepeat(tempRepIndex) = tempRepeat(tempRepIndex) + 1;        % add to that bin in [repeat]
                % confidence analysis on only repeat choices
                tempRepConfTot(tempRepIndex) = tempRepConfTot(tempRepIndex) + 1;    % count in general
                if expe.rslt(iblock).conf(i) == 2                                   % if confidence is high
                    tempRepConf(tempRepIndex) = tempRepConf(tempRepIndex) + 1;      % count in high confidence
                end
            end
        end
        
        % fill in actual holding structures w/ temp arrays by filtering condition
        if strcmpi(blockFilter,'direction')
            if expe.blck(iblock).taskid == 1 %postdiction
                resps1      = vertcat(resps1, tempRevResp);
                revconfs1   = vertcat(revconfs1, tempRevConf);
                expResps1   = vertcat(expResps1, expeRevResp);
                
                reptot1         = vertcat(reptot1, tempRepTot);
                repeat1         = vertcat(repeat1, tempRepeat);
                repconfsCtr1    = vertcat(repconfsCtr1, tempRepConf);
                repconfsTotCtr1 = vertcat(repconfsTotCtr1, tempRepConfTot);
                
                signedconfsCtr1     = vertcat(signedconfsCtr1, tempSignedConf);
                signedconfsTotCtr1  = vertcat(signedconfsTotCtr1, tempSignedConfTot);
                
                % code if looking at interactions
                if expe.blck(iblock).condtn == 1 & interaction% past-low
                    resps1_lo      = vertcat(resps1_lo, tempRevResp);
                    revconfs1_lo   = vertcat(revconfs1_lo, tempRevConf);
                    expResps1_lo   = vertcat(expResps1_lo, expeRevResp);

                    reptot1_lo         = vertcat(reptot1_lo, tempRepTot);
                    repeat1_lo         = vertcat(repeat1_lo, tempRepeat);
                    repconfsCtr1_lo    = vertcat(repconfsCtr1_lo, tempRepConf);
                    repconfsTotCtr1_lo = vertcat(repconfsTotCtr1_lo, tempRepConfTot);

                    signedconfsCtr1_lo     = vertcat(signedconfsCtr1_lo, tempSignedConf);
                    signedconfsTotCtr1_lo  = vertcat(signedconfsTotCtr1_lo, tempSignedConfTot);
                elseif expe.blck(iblock).condtn == 2 & interaction% past-high
                    resps1_hi      = vertcat(resps1_hi, tempRevResp);
                    revconfs1_hi   = vertcat(revconfs1_hi, tempRevConf);
                    expResps1_hi   = vertcat(expResps1_hi, expeRevResp);

                    reptot1_hi         = vertcat(reptot1_hi, tempRepTot);
                    repeat1_hi         = vertcat(repeat1_hi, tempRepeat);
                    repconfsCtr1_hi    = vertcat(repconfsCtr1_hi, tempRepConf);
                    repconfsTotCtr1_hi = vertcat(repconfsTotCtr1_hi, tempRepConfTot);

                    signedconfsCtr1_hi     = vertcat(signedconfsCtr1_hi, tempSignedConf);
                    signedconfsTotCtr1_hi  = vertcat(signedconfsTotCtr1_hi, tempSignedConfTot);
                elseif interaction
                    error(['Unexpected condtn value: ' num2str(expe.blck(iblock).condtn)]);
                end
                
            elseif expe.blck(iblock).taskid == 2 %prediction
                resps2      = vertcat(resps2, tempRevResp);
                revconfs2   = vertcat(revconfs2, tempRevConf);
                expResps2   = vertcat(expResps2, expeRevResp);
                
                reptot2         = vertcat(reptot2, tempRepTot);
                repeat2         = vertcat(repeat2, tempRepeat);
                repconfsCtr2    = vertcat(repconfsCtr2, tempRepConf);
                repconfsTotCtr2 = vertcat(repconfsTotCtr2, tempRepConfTot);
                
                signedconfsCtr2     = vertcat(signedconfsCtr2, tempSignedConf);
                signedconfsTotCtr2  = vertcat(signedconfsTotCtr2, tempSignedConfTot);
                
                % code if looking at interactions
                if expe.blck(iblock).condtn == 1 & interaction% future-low
                    resps2_lo      = vertcat(resps2_lo, tempRevResp);
                    revconfs2_lo   = vertcat(revconfs2_lo, tempRevConf);
                    expResps2_lo   = vertcat(expResps2_lo, expeRevResp);

                    reptot2_lo         = vertcat(reptot2_lo, tempRepTot);
                    repeat2_lo         = vertcat(repeat2_lo, tempRepeat);
                    repconfsCtr2_lo    = vertcat(repconfsCtr2_lo, tempRepConf);
                    repconfsTotCtr2_lo = vertcat(repconfsTotCtr2_lo, tempRepConfTot);

                    signedconfsCtr2_lo     = vertcat(signedconfsCtr2_lo, tempSignedConf);
                    signedconfsTotCtr2_lo  = vertcat(signedconfsTotCtr2_lo, tempSignedConfTot);
                elseif expe.blck(iblock).condtn == 2 & interaction% future-high
                    resps2_hi      = vertcat(resps2_hi, tempRevResp);
                    revconfs2_hi   = vertcat(revconfs2_hi, tempRevConf);
                    expResps2_hi   = vertcat(expResps2_hi, expeRevResp);

                    reptot2_hi         = vertcat(reptot2_hi, tempRepTot);
                    repeat2_hi         = vertcat(repeat2_hi, tempRepeat);
                    repconfsCtr2_hi    = vertcat(repconfsCtr2_hi, tempRepConf);
                    repconfsTotCtr2_hi = vertcat(repconfsTotCtr2_hi, tempRepConfTot);

                    signedconfsCtr2_hi     = vertcat(signedconfsCtr2_hi, tempSignedConf);
                    signedconfsTotCtr2_hi  = vertcat(signedconfsTotCtr2_hi, tempSignedConfTot);
                elseif interaction
                    error(['Unexpected condtn value: ' num2str(expe.blck(iblock).condtn)]);
                end
                
            else
                error(['unexpected value of taskid on block ' num2str(iblock)]);
            end
        elseif strcmpi(blockFilter, 'volatility')
            if expe.blck(iblock).condtn == 1
                resps1      = vertcat(resps1, tempRevResp);
                revconfs1   = vertcat(revconfs1, tempRevConf);
                expResps1   = vertcat(expResps1, expeRevResp);
                
                reptot1         = vertcat(reptot1, tempRepTot);
                repeat1         = vertcat(repeat1, tempRepeat);
                repconfsCtr1    = vertcat(repconfsCtr1, tempRepConf);
                repconfsTotCtr1 = vertcat(repconfsTotCtr1, tempRepConfTot);
                
                signedconfsCtr1     = vertcat(signedconfsCtr1, tempSignedConf);
                signedconfsTotCtr1  = vertcat(signedconfsTotCtr1, tempSignedConfTot);
            elseif expe.blck(iblock).condtn == 2
                resps2      = vertcat(resps2, tempRevResp);
                revconfs2   = vertcat(revconfs2, tempRevConf);
                expResps2   = vertcat(expResps2, expeRevResp);
                
                reptot2         = vertcat(reptot2, tempRepTot);
                repeat2         = vertcat(repeat2, tempRepeat);
                repconfsCtr2    = vertcat(repconfsCtr2, tempRepConf);
                repconfsTotCtr2 = vertcat(repconfsTotCtr2, tempRepConfTot);
                
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
        subjIndexTracker = vertcat(subjIndexTracker, size(resps1,1)+1);
    end
end
clear expe expeResp i iblock isubj subject subjChoice subjResp filename;

%% Reversal Curves
clf;
close all;
% calculate rate of responses (correct after the change) around change point
diff1 = abs(expResps1-resps1);
diff2 = abs(expResps2-resps2);
nChanges1 = numel(expResps1(:,1));
nChanges2 = numel(expResps2(:,1));
% rates for Figure 11
rates1_1 = zeros(1,befAftTrials*2);
rates1_2 = zeros(1,befAftTrials*2);
confrates1_1 = zeros(1,befAftTrials*2);
confrates1_2 = zeros(1,befAftTrials*2);

% rates on individual subjects
ssrate1 = zeros(1,8);
ssrate2 = zeros(1,8);

ssconf1 = zeros(1,8);
ssconf2 = zeros(1,8);

subjCtr = 1;
for i = 1:numel(resps1(:,1))
    ssrate1(subjCtr,:) = ssrate1(subjCtr,:) + diff1(i,:);
    ssrate2(subjCtr,:) = ssrate2(subjCtr,:) + diff2(i,:);
    
    modNum = numel(diff1(:,1))./25;
        
    if mod(i,modNum) == 0
        ssrate1(subjCtr,:) = ssrate1(subjCtr,:)./modNum;
        ssrate2(subjCtr,:) = ssrate2(subjCtr,:)./modNum;

        ssconf1(subjCtr,:) = sum(revconfs1(i-modNum+1:i,:)-ones(modNum,8))./modNum;
        ssconf2(subjCtr,:) = sum(revconfs2(i-modNum+1:i,:)-ones(modNum,8))./modNum;
        
        if i ~= numel(resps1(:,1))
            ssrate1 = vertcat(ssrate1,zeros(1,8));
            ssrate2 = vertcat(ssrate2,zeros(1,8));
        end
        subjCtr = subjCtr + 1;
    end
end

% statistics on repetition rate in the two conditions
[hRev statsRev] = ttest(ssrate1, ssrate2);
statsRevAsterisks = {};
statCtr = 1;
for stat = statsRev
    if stat < .001 
        statsRevAsterisks = [statsRevAsterisks '***'];
    elseif stat < .01
        statsRevAsterisks = [statsRevAsterisks '**'];
    elseif stat < .05
        statsRevAsterisks = [statsRevAsterisks '*'];
    else
        statsRevAsterisks = [statsRevAsterisks ' '];
    end
    statCtr = statCtr + 1;
end

% rates across all subjects
for i = 1:befAftTrials*2
    rates1_1(i)     = numel(diff1(diff1(:,i)==1))/nChanges1;
    rates1_2(i)     = numel(diff2(diff2(:,i)==1))/nChanges2;
    confrates1_1(i) = sum(revconfs1(:,i)==2)/nChanges1;
    confrates1_2(i) = sum(revconfs2(:,i)==2)/nChanges2;
end

% getting stuff from the model simulations
load('modelStruct');
rev_det = [];
rev_inf = [];
rev_sel = [];
revconf_det = [];
revconf_inf = [];
revconf_sel = [];
for i = 1:30
    if ismember(i, excluded)
       continue; 
    end
    rev_det = vertcat(rev_det, modelStruct(i).det.past.rev);
    rev_inf = vertcat(rev_inf, modelStruct(i).inf.past.rev);
    rev_sel = vertcat(rev_sel, modelStruct(i).sel.past.rev);
    revconf_det = vertcat(revconf_det, modelStruct(i).det.past.revconf);
    revconf_inf = vertcat(revconf_inf, modelStruct(i).inf.past.revconf);
    revconf_sel = vertcat(revconf_sel, modelStruct(i).sel.past.revconf);
end

% interaction stuff
if interaction
    % convert subj responses for the learning rate calc.
    diff1_lo = abs(expResps1_lo-resps1_lo);
    diff1_hi = abs(expResps1_hi-resps1_hi);
    diff2_lo = abs(expResps2_lo-resps2_lo);
    diff2_hi = abs(expResps2_hi-resps2_hi);
    
    rev1_lo = [];
    rev1_hi = [];
    rev2_lo = [];
    rev2_hi = [];
    
    modNum_lo = numel(diff1_lo(:,1))./25;
    modNum_hi = numel(diff1_hi(:,1))./25;
    % loop through low volatility stuff
    for i = 1:numel(resps1_lo(:,1))
        if mod(i,modNum_lo) == 0
            rev1_lo = vertcat(rev1_lo, sum(diff1_lo(i-modNum_lo+1:i,:))./modNum_lo);
            rev2_lo = vertcat(rev2_lo, sum(diff2_lo(i-modNum_lo+1:i,:))./modNum_lo);
        end
    end
    % loop through high volatility stuff
    for i = 1:numel(resps1_hi(:,1))
        if mod(i,modNum_hi) == 0
            rev1_hi = vertcat(rev1_hi, sum(diff1_hi(i-modNum_hi+1:i,:))./modNum_hi);
            rev2_hi = vertcat(rev2_hi, sum(diff2_hi(i-modNum_hi+1:i,:))./modNum_hi);
        end
    end   
end

%% Plots: Interaction : Inferential Direction x Volatility

if interaction
    
    ttest(rev1_lo, rev1_hi)
    ttest(rev2_lo, rev2_hi)
    
    axisSupp1 = [-befAftTrials:befAftTrials];
    axisSupp1 = axisSupp1(axisSupp1~=0);  
    figure(1);
    hold on;
    title(['Postdiction x Low/High Volatility']);
    p1_1 = scatter(axisSupp1, mean(rev1_lo),'MarkerEdgeColor',[0, 0.4470, 0.7410]);
    errorbar(axisSupp1, mean(rev1_lo),std(rev1_lo)./sqrt(25),'LineStyle','none','Color',[0, 0.4470, 0.7410]);
    p1_2 = scatter(axisSupp1, mean(rev1_hi),'MarkerEdgeColor',[0.8500, 0.3250, 0.0980]);
    errorbar(axisSupp1, mean(rev1_hi),std(rev1_hi)./sqrt(25),'LineStyle','none','Color',[0.8500, 0.3250, 0.0980]);
    legend([p1_1 p1_2],{'Low','High'},'Location','southeast');
    hold off;
    figure(2);
    hold on;
    title(['Prediction x Low/High Volatility']);
    p2_1 = scatter(axisSupp1, mean(rev2_lo),'MarkerEdgeColor',[0, 0.4470, 0.7410]);
    errorbar(axisSupp1, mean(rev2_lo),std(rev2_hi)./sqrt(25),'LineStyle','none','Color',[0, 0.4470, 0.7410]);
    p2_2 = scatter(axisSupp1, mean(rev2_hi),'MarkerEdgeColor',[0.8500, 0.3250, 0.0980]);
    errorbar(axisSupp1, mean(rev2_hi),std(rev2_hi)./sqrt(25),'LineStyle','none','Color',[0.8500, 0.3250, 0.0980]);
    legend([p2_1 p2_2],{'Low','High'},'Location','southeast');
    hold off; 
end


%% Plots: Reversal Curves : Learning 

data2modelComp = true; % set to true to only check condition 1 data to condition 1 models

axisSupp1 = [-befAftTrials:befAftTrials];
axisSupp1 = axisSupp1(axisSupp1~=0);
% This plots the percentage of responses corresponding to the CORRECT STATE of the
% environment AFTER the switch point
figure(10);
hold on;
%models
if data2modelComp
    %stats
    [p h] = ttest(ssrate1, rev_det);
    h = [h(1:befAftTrials) NaN h(befAftTrials+1:end)];
    asts = asterisks(h);
    plotAsterisks([-4:4], .95, asts, detRGB);
    [p h] = ttest(ssrate1, rev_inf);
    h = [h(1:befAftTrials) NaN h(befAftTrials+1:end)];
    asts = asterisks(h);
    plotAsterisks([-4:4], .93, asts, infRGB);
    [p h] = ttest(ssrate1, rev_sel);
    h = [h(1:befAftTrials) NaN h(befAftTrials+1:end)];
    asts = asterisks(h);
    plotAsterisks([-4:4], .91, asts, selRGB);
    
    p10_1det        = plot(axisSupp1,mean(rev_det),'Color',detRGB);
    p10_1det_err    = shadedErrorBar(axisSupp1, mean(rev_det), std(rev_det)./sqrt(25),'lineprops',{'--','Color',detRGB},'transparent',1);
    p10_1inf        = plot(axisSupp1,mean(rev_inf),'Color',infRGB);
    p10_1inf_err    = shadedErrorBar(axisSupp1, mean(rev_inf), std(rev_inf)./sqrt(25),'lineprops',{'--','Color',infRGB},'transparent',1);
    p10_1sel        = plot(axisSupp1,mean(rev_sel),'Color',selRGB);
    p10_1sel_err    = shadedErrorBar(axisSupp1, mean(rev_sel), std(rev_sel)./sqrt(25),'lineprops',{'--','Color',selRGB},'transparent',1);
end
%data
p10_1a = scatter(axisSupp1, mean(ssrate1),'MarkerEdgeColor',postRGB,'MarkerFaceColor',postRGB);
p10_1s = errorbar(axisSupp1, mean(ssrate1),std(ssrate1)/sqrt(length(ssrate1)),'LineStyle','none',...
                'Color',postRGB);    % condition 1 stdev over single subject means
if ~data2modelComp
    p10_2a = scatter(axisSupp1, mean(ssrate2),'MarkerEdgeColor',predRGB,'MarkerFaceColor',predRGB);
    p10_2s = errorbar(axisSupp1, mean(ssrate2),std(ssrate2)/sqrt(length(ssrate2)),'LineStyle','none',...
                    'Color',predRGB);   % condition 2 stdev over single subject means
end
xline(0);
xlim([-4.2 4.2]);
xlabel('Trial position before/after reversal point (i.e. 0)');
ylim([0 1]);
if data2modelComp
    legend([p10_1a p10_1det p10_1inf p10_1sel], {txtCond1,'Model: det','Model: inf','Model: inf+sel'},'Location','southeast');
    title(sprintf('Comparison of behavior (models to data) on reversal learning in the %s condition',blockFilter));
else
    asterisksY = mean(ssrate2);     % statistical significance over interval
    for i = 1:befAftTrials*2
        if i<5
            i2 = i-5;
        else 
            i2 = i-4;
        end
        text(i2-.05, asterisksY(i)+.08, cell2mat(statsRevAsterisks(i)),'FontSize',14);
    end  
    legend([p10_1a p10_2a], {txtCond1, txtCond2},'Location','southeast');
    title({sprintf('Effect of %s on reversal learning',blockFilter) '25 subjects; error bars SEM'});
end
hold off;

%% Plots: Reversal Curves : Confidence

figure(11);
hold on;
%models
if data2modelComp
    %stats
    [p h] = ttest(ssconf1, revconf_det);
    h = [h(1:befAftTrials) NaN h(befAftTrials+1:end)];
    asts = asterisks(h);
    plotAsterisks([-4:4], .64, asts, detRGB);
    
    [p h] = ttest(ssconf1, revconf_inf);
    h = [h(1:befAftTrials) NaN h(befAftTrials+1:end)];
    asts = asterisks(h);
    plotAsterisks([-4:4], .62, asts, infRGB);
    
    [p h] = ttest(ssconf1, revconf_sel);
    h = [h(1:befAftTrials) NaN h(befAftTrials+1:end)];
    asts = asterisks(h);
    plotAsterisks([-4:4], .6, asts, selRGB);
    
    p11_1detconf        = plot(axisSupp1,mean(revconf_det),'Color',detRGB);
    p11_1detconf_err    = shadedErrorBar(axisSupp1, mean(revconf_det), std(revconf_det)./sqrt(25),'lineprops',{'--','Color',detRGB},'transparent',1,...
                                         'patchSaturation',.1);
    p11_1infconf        = plot(axisSupp1,mean(revconf_inf),'Color',infRGB);
    p11_1infconf_err    = shadedErrorBar(axisSupp1, mean(revconf_inf), std(revconf_inf)./sqrt(25),'lineprops',{'--','Color',infRGB},'transparent',1,...
                                         'patchSaturation',.1);
    p11_1selconf        = plot(axisSupp1,mean(revconf_sel),'Color',selRGB);
    p11_1selconf_err    = shadedErrorBar(axisSupp1, mean(revconf_sel), std(revconf_sel)./sqrt(25),'lineprops',{'--','Color',selRGB},'transparent',1,...
                                         'patchSaturation',.1);
end
%data
p11_1s = errorbar(axisSupp1, mean(ssconf1),std(ssconf1)/sqrt(length(ssconf1)),'LineStyle','none',...
                'Color',postRGB);    % condition 1 stdev over single subject means
p11_1a = scatter(axisSupp1,mean(ssconf1),'MarkerEdgeColor',postRGB,'MarkerFaceColor',postRGB);

if ~data2modelComp
    p11_2s = errorbar(axisSupp1, mean(ssconf2),std(ssconf2)/sqrt(length(ssconf2)),'LineStyle','none',...
                    'Color',predRGB);    % condition 2 stdev over single subject means
    p11_2a = scatter(axisSupp1,mean(ssconf2),'MarkerEdgeColor',predRGB,'MarkerFaceColor',predRGB);
end

xline(0);
yline(0.5,':');
xlim([-4.1 4.1]);
xlabel('Trial position before/after reversal point (i.e. 0)');
if data2modelComp
    legend([p11_1a p11_1detconf p11_1infconf p11_1selconf], {txtCond1,'Model: det','Model: inf','Model: inf+sel'},'Location','southeast');
    title(sprintf('Comparison of behavior (models to data) on rates of high confidence \n during reversals in the %s condition',blockFilter));
else
    title({sprintf('Effect of %s on confidence around reversal point',blockFilter) '25 subjects; error bars SEM'});
    legend([p11_1a p11_2a], {txtCond1, txtCond2},'Location','southeast');
end
hold off;


%% Percentage of repeated choices over seqLLR signed by choice
% calculate rate of repeated choices from counters in main loop

% means
reprate1 = zeros(1,8);
reprate2 = zeros(1,8);
repconfrate1 = zeros(1,8);
repconfrate2 = zeros(1,8);
signedconfrate1 = zeros(1,8);
signedconfrate2 = zeros(1,8);

for i = 1:8
    % mean calculations
    reprate1(i) = sum(repeat1(:,i))/sum(reptot1(:,i));
    reprate2(i) = sum(repeat2(:,i))/sum(reptot2(:,i));
    repconfrate1(i) = sum(repconfsCtr1(:,i))/sum(repconfsTotCtr1(:,i));
    repconfrate2(i) = sum(repconfsCtr2(:,i))/sum(repconfsTotCtr2(:,i));
    signedconfrate1(i) = sum(signedconfsCtr1(:,i))/sum(signedconfsTotCtr1(:,i));
    signedconfrate2(i) = sum(signedconfsCtr2(:,i))/sum(signedconfsTotCtr2(:,i));
end

% The following code is not optimized for text space. It is separated into redundant
% loops for clarity.

% data organized to generate statistics for subject-wide variability on REPETITION RATE
sssignreprate1 = zeros(1,8);
sssignreprate2 = zeros(1,8);
subjCtr = 1;
for i = 1:numel(repeat1(:,1))
        sssignreprate1(subjCtr,:) = sssignreprate1(subjCtr,:) + repeat1(i,:);
        sssignreprate2(subjCtr,:) = sssignreprate2(subjCtr,:) + repeat2(i,:);
    if mod(i,4) == 0
        sssignreprate1(subjCtr,:) = sssignreprate1(subjCtr,:)./sum(reptot1(i-3:i,:));
        sssignreprate2(subjCtr,:) = sssignreprate2(subjCtr,:)./sum(reptot2(i-3:i,:));
        if i ~= numel(repeat1(:,1))
            sssignreprate1 = vertcat(sssignreprate1,zeros(1,8));
            sssignreprate2 = vertcat(sssignreprate2,zeros(1,8));
        end
        subjCtr = subjCtr + 1;
    end
end

% statistics on repetition rate in the two conditions
[hRepSigned statsRepSigned] = ttest(sssignreprate1, sssignreprate2);
statsRepSignedAsterisks = {};
statCtr = 1;
for stat = statsRepSigned
    if stat < .001 
        statsRepSignedAsterisks = [statsRepSignedAsterisks '***'];
    elseif stat < .01
        statsRepSignedAsterisks = [statsRepSignedAsterisks '**'];
    elseif stat < .05
        statsRepSignedAsterisks = [statsRepSignedAsterisks '*'];
    else
        statsRepSignedAsterisks = [statsRepSignedAsterisks ' '];
    end
    statCtr = statCtr + 1;
end

% data organized to generate statistics for subject-wide variability on CONFIDENCE
sssignconfrate1 = zeros(1,8); % single subject rate of high confidence on condition 1
sssignconfrate2 = zeros(1,8); % single subject rate of high confidence on condition 2
subjCtr = 1;
for i = 1:numel(signedconfsCtr1(:,1))
        sssignconfrate1(subjCtr,:) = sssignconfrate1(subjCtr,:) + signedconfsCtr1(i,:);
        sssignconfrate2(subjCtr,:) = sssignconfrate2(subjCtr,:) + signedconfsCtr2(i,:);
    if mod(i,4) == 0
        sssignconfrate1(subjCtr,:) = sssignconfrate1(subjCtr,:)./sum(reptot1(i-3:i,:));
        sssignconfrate2(subjCtr,:) = sssignconfrate2(subjCtr,:)./sum(reptot2(i-3:i,:));
        if i ~= numel(signedconfsCtr1(:,1))
            sssignconfrate1 = vertcat(sssignconfrate1,zeros(1,8));
            sssignconfrate2 = vertcat(sssignconfrate2,zeros(1,8));
        end
        subjCtr = subjCtr + 1;
    end
end

% statistics on confidence in the two conditions
[hConfSigned statsConfSigned] = ttest(sssignconfrate1, sssignconfrate2);
statsConfSignedAsterisks = {};
statCtr = 1;
for stat = statsConfSigned
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

% getting stuff from the model simulations
load('modelStruct');
rep_det = [];
rep_inf = [];
rep_sel = [];
signconf_det = [];
signconf_inf = [];
signconf_sel = [];
%model
for i = 1:30
    if ismember(i, excluded)
       continue; 
    end
    rep_det = vertcat(rep_det, modelStruct(i).det.past.rep);
    rep_inf = vertcat(rep_inf, modelStruct(i).inf.past.rep);
    rep_sel = vertcat(rep_sel, modelStruct(i).sel.past.rep);
    signconf_det = vertcat(signconf_det, modelStruct(i).det.past.signconf);
    signconf_inf = vertcat(signconf_inf, modelStruct(i).inf.past.signconf);
    signconf_sel = vertcat(signconf_sel, modelStruct(i).sel.past.signconf);
end


%% This plots the percentage of repeated choices over signed seqllr
figure(20);
hold on;
%models  
if data2modelComp
    %stats
    [p h] = ttest(sssignreprate1, rep_det)
    h = [h(1:befAftTrials) NaN h(befAftTrials+1:end)];
    asts = asterisks(h);
    plotAsterisks([1:8], .64, asts, detRGB);
    
    [p h] = ttest(sssignreprate1, rep_inf)
    h = [h(1:befAftTrials) NaN h(befAftTrials+1:end)];
    asts = asterisks(h);
    plotAsterisks([1:8], .62, asts, infRGB);
    
    [p h] = ttest(sssignreprate1, rep_sel)
    h = [h(1:befAftTrials) NaN h(befAftTrials+1:end)];
    asts = asterisks(h);
    plotAsterisks([1:8], .6, asts, selRGB);
    
    p20_1det        = plot([1:8],mean(rep_det),'Color',detRGB);
    p20_1det_err    = shadedErrorBar([1:8], mean(rep_det), std(rep_det)./sqrt(25),'lineprops',{'--','Color',detRGB},'transparent',1);
    p20_1inf        = plot([1:8],mean(rep_inf),'Color',infRGB);
    p20_1inf_err    = shadedErrorBar([1:8], mean(rep_inf), std(rep_inf)./sqrt(25),'lineprops',{'--','Color',infRGB},'transparent',1);
    p20_1sel        = plot([1:8],mean(rep_sel),'Color',selRGB);
    p20_1inf_err    = shadedErrorBar([1:8], mean(rep_sel), std(rep_sel)./sqrt(25),'lineprops',{'--','Color',selRGB},'transparent',1);
end
%data
% condition 1 mean over single subject means
p20_1a = scatter([1:8], mean(sssignreprate1),'MarkerEdgeColor',postRGB,'LineWidth',2','MarkerFaceColor',postRGB);    
p20_1s = errorbar([1:8], mean(sssignreprate1),std(sssignreprate1)/sqrt(length(sssignreprate1)),'LineStyle','none',...
                'Color',postRGB);                                           % condition 1 stdev over single subject means
if ~data2modelComp
    p20_2a = scatter([1:8], mean(sssignreprate2),'MarkerEdgeColor',predRGB);% condition 2 mean over single subject means
    p20_2s = errorbar([1:8], mean(sssignreprate2),std(sssignreprate2)/sqrt(length(sssignreprate2)),'LineStyle','none',...
                    'Color',predRGB);                                       % condition 2 stdev over single subject means
end
     
xticks([1:8]);
xticklabels({'(-\infty,-3)', '[-3:-2)', '[-2:-1)', '[-1:0)', '[0:1)', '[1:2)', '[2:3)', '[3:\infty)'}); 
yline(.5,':');
ylim([0 1]);
if data2modelComp
    legend([p20_1a p20_1det p20_1inf p20_1sel], {txtCond1,'Model: det','Model: inf','Model: inf+sel'},'Location','southeast');
    title(sprintf('Comparison of behavior (models to data) on repeat choices over \n signed evidence in the %s condition',blockFilter));
else
    asterisksY = mean(sssignreprate1);     % statistical significance over interval
    for i = 1:numel(statsRepSignedAsterisks)
        text(i, asterisksY(i)-.08, cell2mat(statsRepSignedAsterisks(i)),'FontSize',14);
    end   
    legend([p20_1a p20_2a], {txtCond1, txtCond2},'Location','southeast');
    title({sprintf('Effect of %s on repeat choices against evidence signed by the previous choice', blockFilter) '25 subjects; error bars SEM'});
end
xlabel('Sequence LLR at trial signed by previous choice');
hold off;

%% **** This code needs rework****
% This plots the percentage of high confidence over signed seqllr on REPEATED choice trials
%{
figure(22);
p32 = plot([1:8], repconfrate1, 'LineWidth', 2);
hold on;
p42 = plot([1:8], repconfrate2, 'LineWidth', 2);
title(sprintf('Effect of %s on confidence on repeat choice trials vs evidence signed by previous choice', blockFilter));
xticks([1:8]);
xticklabels({'(-\infty,-3)', '[-3:-2)', '[-2:-1)', '[-1:0)', '[0:1)', '[1:2)', '[2:3)', '[3:\infty)'}); 
legend([p32,p42], {txtCond1, txtCond2},'Location','southeast');
hold off;
%}

%% This plots the percentage of high confidence over signed seqllr on ALL TRIALS
figure(30);
hold on;
%model
if data2modelComp
    %stats
    [p h] = ttest(sssignconfrate1, signconf_det)
    h = [h(1:befAftTrials) NaN h(befAftTrials+1:end)];
    asts = asterisks(h);
    plotAsterisks([1:8], .95, asts, detRGB);
    
    [p h] = ttest(sssignconfrate1, signconf_inf)
    h = [h(1:befAftTrials) NaN h(befAftTrials+1:end)];
    asts = asterisks(h);
    plotAsterisks([1:8], .93, asts, infRGB);
    
    [p h] = ttest(sssignconfrate1, signconf_sel)
    h = [h(1:befAftTrials) NaN h(befAftTrials+1:end)];
    asts = asterisks(h);
    plotAsterisks([1:8], .91, asts, selRGB);
    p30_1detsignconf       = plot([1:8],mean(signconf_det),'Color',detRGB);
    p30_1detsignconf_err   = shadedErrorBar([1:8], mean(signconf_det), std(signconf_det)./sqrt(25),'lineprops',{'--','Color',detRGB},'transparent',1,...
                                        'patchSaturation',.1);
    p30_1infsignconf       = plot([1:8],mean(signconf_inf),'Color',infRGB);
    p30_1infsignconf_err   = shadedErrorBar([1:8], mean(signconf_inf), std(signconf_inf)./sqrt(25),'lineprops',{'--','Color',infRGB},'transparent',1,...
                                         'patchSaturation',.1);
    p30_1selsignconf       = plot([1:8],mean(signconf_sel),'Color',selRGB);
    p30_1selsignconf_err   = shadedErrorBar([1:8], mean(signconf_sel), std(signconf_sel)./sqrt(25),'lineprops',{'--','Color',selRGB},'transparent',1,...
                                         'patchSaturation',.1);
end

%data
p30_1s = errorbar([1:8], mean(sssignconfrate1),std(sssignconfrate1)/sqrt(length(sssignconfrate1)),'LineStyle','none',...
                'Color',[1 .75 .5]);    % condition 1 stdev over single subject means
p30_1a = scatter([1:8], mean(sssignconfrate1),'MarkerEdgeColor',postRGB);  % condition 1 mean over single subject means
if ~data2modelComp
    p30_1s = errorbar([1:8], mean(sssignconfrate2),std(sssignconfrate2)/sqrt(length(sssignconfrate2)),'LineStyle','none',...
                    'Color',[.75 .75 1]);   % condition 2 stdev over single subject meansj
    p30_1a = scatter([1:8], mean(sssignconfrate2),'MarkerEdgeColor',predRGB); % condition 2 mean over single subject means
end

xlim([.9 8.1]);
xticks([1:8]);
xticklabels({'(-\infty,-3)', '[-3:-2)', '[-2:-1)', '[-1:0)', '[0:1)', '[1:2)', '[2:3)', '[3:\infty)'}); 
xlabel('Sequence LLR at trial $\textit{t}$ signed by choice at trial $\textit{t-1}$','Interpreter','latex');

ylim([0 1]);
if ~data2modelComp
    asterisksY = mean(sssignconfrate1);     % statistical significance over interval
    for i = 1:numel(statsConfSignedAsterisks)
        text(i, asterisksY(i)+.08, cell2mat(statsConfSignedAsterisks(i)),'FontSize',14);
    end
    title({sprintf('Effect of %s on confidence vs evidence signed by the previous choice', blockFilter) '25 subjects; error bars SEM'});
    legend([p30_1a,p30_2], {txtCond1, txtCond2},'Location','southeast');
else
    title(sprintf('Comparison of behavior (models to data) on confidence during on choices over \n signed evidence in the %s condition',blockFilter));
    legend([p30_1a,p30_1detsignconf,p30_1infsignconf,p30_1selsignconf], {txtCond1,'Model: det','Model: inf','Model: inf+sel'},'Location','southeast');
    
end
hold off;

%% Analysis of Unexcluded vs. Excluded Subjects Comparison Code
%   Read carefully before trying to run this section of code
%   1. Run the setup and all other plots on ALL subjects FIRST and do not run this
%       section, or run the whole thing and it will throw an error
%   2. Run the setup and all other plots on subjects after exclusion
%       and then run this section 
if skipExclComp 
    disp('Skipping comparison between subjects that are excluded and not...');
else
    if exist('revCurvePost')==1 && exist('revCurvePost_excl')==1
        figure(3);
        cp13 = plot(repCurveXAxis, repCurvePost);
        hold on;
        cp23 = plot(repCurveXAxis, repCurvePost_excl);
        title('Postdiction repetition curves');
        xticks([1:8]);
        xticklabels({'(-\infty,-3)', '[-3:-2)', '[-2:-1)', '[-1:0)', '[0:1)', '[1:2)', '[2:3)', '[3:\infty)'}); 
        legend([cp13,cp23], {'all', 'with exclusions'},'Location','southeast');
        hold off;

        figure(4);
        cp14 = plot(repCurveXAxis, repCurvePred);
        hold on;
        cp24 = plot(repCurveXAxis, repCurvePred_excl);
        title('Prediction repetition curves');
        xticks([1:8]);
        xticklabels({'(-\infty,-3)', '[-3:-2)', '[-2:-1)', '[-1:0)', '[0:1)', '[1:2)', '[2:3)', '[3:\infty)'}); 
        legend([cp14,cp24], {'all', 'with exclusions'},'Location','southeast');
        hold off;

        figure(5);
        cp15 = plot(revCurveXAxis, revCurvePred);
        hold on;
        cp25 = plot(revCurveXAxis, revCurvePred_excl);
        title('Prediction reversal curves');
        legend([cp15,cp25], {'all', 'with exclusions'},'Location','southeast');
        hold off;

        figure(6);
        cp16 = plot(revCurveXAxis, revCurvePost);
        hold on;
        cp26 = plot(revCurveXAxis, revCurvePost_excl);
        title('Postdiction reversal curves');
        legend([cp16,cp26], {'all', 'with exclusions'},'Location','southeast');
        hold off;
    else
        error('Cannot find all relevant data. Make sure that you''ve run the analysis for both exclusion conditions!');
    end
end

clear cp13 cp14 cp15 cp16 cp23 cp24 cp25 cp26;

%% Local Functions

function astArray = asterisks(hvals)
    astArray = {};
    for hval = hvals
        if hval < .001 
            astArray = [astArray '***'];
        elseif hval < .01
            astArray = [astArray '**'];
        elseif hval < .05
            astArray = [astArray '*'];
        else
            astArray = [astArray ' '];
        end
    end
end

function plotAsterisks(xsupport, y, asterisks, color)
    for i = 1:numel(xsupport)
        ast = cell2mat(asterisks(i));
        text(xsupport(i),y,ast,'Color', color);
    end
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
