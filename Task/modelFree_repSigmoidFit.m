% Model free analysis on TURFU subjects
%   (Sigmoid curve fitting on individual subject repeated choice curves on signed
%    seqllr)
%
%   Name:   modelFree_repSigmoidFit.m
%   Type:   Script
% Author:   Jun Seok Lee
%   Date:   April 2019
% 
% Description:  Creates curves of repeat choice percentage vs evidence signed by previous choice
%               for each subject.
%               Fits a sigmoid curve for each subject on the two conditions.

clf;
close all;
% a whole lotta colors
postRGB = [1 .75 .5];
predRGB = [.75 .75 1];
%% Setup
% Static variables (set and forget)
n_subjects  = 30;   % set to number of subjects we are analyzing
excluded    = [21]; % excluded subject numbers

% Dynamic variables (CHANGE VALUES HERE BASED ON ANALYSIS)
blockFilter     = 'direction'; % 'direction' or 'volatility'
exclude         = true;        % set to true to exclude shitty results
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

% Subject responses and confidence values (organized by filter)
respTot1 = zeros(n_subjects,8);
respTot2 = zeros(n_subjects,8);
respRepeat1 = zeros(n_subjects,8);
respRepeat2 = zeros(n_subjects,8);
repeats1 = []; % rates of repetitive choice given signed seqllr (8 bins) in 1st condition
repeats2 = []; % rates of repetitive choice given signed seqllr (8 bins) in 2nd condition

% Parameters for the logistic regression fit on the two conditions
postParams = zeros(n_subjects,2);    % sigmoid curve parameters on the postdictive condition
predParams = zeros(n_subjects,2);    % sigmoid curve parameters on the predictive condition

%% Data point (rate of repeptive choice vs signed seqllr) production for each subject

% General algorithm:
%
%   FOR each subject
%       FOR all real blocks
%           get relevant data 
%           organize by condition
%
%       get data points (rates) for each subject with the appropriate bin
%
%       FOR each condition
%           fit a sigmoid curve that fits the condition
%           store the parameters of curve in holding data structure

for isubj = 1:n_subjects
    
    if ismember(isubj, excluded)
        continue;
    end
    
    %load subject information
    subject  = isubj;
    filename = dir(sprintf('Data/TURFU_S%02d_*.mat', subject));
    filename = filename.name;
    load(sprintf('Data/%s',filename));
    
    % set up subject-temporary structures
    tempSubjRespsTot1    = zeros(1,8);
    tempSubjRespsRepeat1 = zeros(1,8);
    tempSubjRespsTot2    = zeros(1,8);
    tempSubjRespsRepeat2 = zeros(1,8);
    
    % go through the blocks
    for iblock = 3:10
        % set up block-temporary structures to be filled into higher level subj structs
        % ## 1 to 8 corresponds to bin indices described in the function (repeatBin) below 
        tempRespsTot    = zeros(1,8);
        tempRespsRepeat = zeros(1,8);
        
        % go through the trials
        for itrial = 3:73 % in the counting scheme for the subject (1 is nothing, 2 is first trial)
            % identify the sign of the subject choice on the previous trial
            prevSubjChoice  = (expe.rslt(iblock).resp(itrial-1)-1.5)*-2; % converts/maps choices (1,2) to (+1,-1)
            tempBinIndex    = repeatBin(prevSubjChoice*expe.blck(iblock).seqllr(itrial-1));
            
            tempRespsTot(tempBinIndex) = tempRespsTot(tempBinIndex) + 1;
            
            currentResp     = expe.rslt(iblock).resp(itrial);
            previousResp    = expe.rslt(iblock).resp(itrial-1);
            if currentResp == previousResp
                tempRespsRepeat(tempBinIndex) = tempRespsRepeat(tempBinIndex) + 1;
            end
        end
        
        % once all trials analyzed in a single block
        % insert block-temp data into subject-temp structure based on filtering condition
        taskid = expe.blck(iblock).taskid;
        if strcmpi(blockFilter, 'direction')
            if taskid == 1
                tempSubjRespsTot1       = tempSubjRespsTot1    + tempRespsTot;
                tempSubjRespsRepeat1    = tempSubjRespsRepeat1 + tempRespsRepeat;
            elseif taskid == 2
                tempSubjRespsTot2       = tempSubjRespsTot2    + tempRespsTot;
                tempSubjRespsRepeat2    = tempSubjRespsRepeat2 + tempRespsRepeat;
            else
                error(['Unexpected value of experimental task id. taskid = ' num2str(taskid)]);
            end
        elseif strcmpi(blockFilter, 'volatility')
            disp('Code for the value ''volatility'' does not exist, yet.');
        else
            error('Invalid value for blockFilter. Set it to ''direction'' or ''volatility''.');
        end
    end
    respTot1(isubj,:) = tempSubjRespsTot1;
    respTot2(isubj,:) = tempSubjRespsTot2;
    respRepeat1(isubj,:) = tempSubjRespsRepeat1;
    respRepeat2(isubj,:) = tempSubjRespsRepeat2;
    
    repeats1(isubj,:) = tempSubjRespsRepeat1./tempSubjRespsTot1;
    repeats2(isubj,:) = tempSubjRespsRepeat2./tempSubjRespsTot2;
    
    % fit logistic curve for each subject
    postParams(isubj,:) = glmfit([1:8], [respRepeat1(isubj,:).' respTot1(isubj,:).'],'binomial','link','logit');
    predParams(isubj,:) = glmfit([1:8], [respRepeat2(isubj,:).' respTot2(isubj,:).'],'binomial','link','logit');
    
end

%get rid of the 0's of excluded subjects
ictr = 1;
for i = 1:30
    if repeats1(ictr,:) == zeros(1,8);
        repeats1(ictr,:) = [];
        repeats2(ictr,:) = [];
        postParams(ictr,:) = [];
        predParams(ictr,:) = [];
        ictr = ictr - 1;
    end
    ictr = ictr + 1;
end

clear previousResp currentResp cond isubj iblock itrial taskid filename expe subject
clearvars -regexp ^temp

%% Plot manipulated data and fitted sigmoid
xAxis = [1:.1:8];
subtractFactor = 0;
figure(1);
hold on;
for isubj = 1:30
    if ismember(isubj, excluded)
        subtractFactor = subtractFactor+1;
        continue;
    end
    isubj = isubj - subtractFactor;

    % logistic curve
    z1(isubj,:) = postParams(isubj,1)+(postParams(isubj,2)*xAxis);
    z1(isubj,:) = 1 ./ (1 + exp(-z1(isubj,:)));
    z2(isubj,:) = predParams(isubj,1)+(predParams(isubj,2)*xAxis);
    z2(isubj,:) = 1 ./ (1 + exp(-z2(isubj,:)));
end
f1_err = shadedErrorBar(xAxis, mean(z1), std(z1)./sqrt(size(z1,1)),'lineprops',{'--','Color',postRGB},'transparent',1,...
               'patchSaturation',.1);
f2_err = shadedErrorBar(xAxis, mean(z2), std(z2)./sqrt(size(z2,1)),'lineprops',{'--','Color',predRGB},'transparent',1,...
               'patchSaturation',.1);
f1 = plot(xAxis,mean(z1),'Color',postRGB);
f2 = plot(xAxis,mean(z2),'Color',predRGB);
errorbar([1:8], mean(repeats1), std(repeats1)/sqrt(size(repeats1,1)),'LineStyle','none','Color',postRGB);
errorbar([1:8], mean(repeats2), std(repeats2)/sqrt(size(repeats2,1)),'LineStyle','none','Color',predRGB);
p1 = scatter([1:8], mean(repeats1),'filled','MarkerEdgeColor',postRGB,'MarkerFaceColor',postRGB);
p2 = scatter([1:8], mean(repeats2),'filled','d','MarkerEdgeColor',predRGB,'MarkerFaceColor',predRGB);
xticks([1:8]);
xticklabels({'(-\infty,-3)', '[-3:-2)', '[-2:-1)', '[-1:0)', '[0:1)', '[1:2)', '[2:3)', '[3:\infty)'}); 
xlabel('Sequence LLR at trial signed by previous choice');

yline(.5,'--');
ylabel('Proportion of repeat choices');
txtSigFit1 = 'Individual sigmoid fits to subject postdiction data';
txtSigFit2 = 'Individual sigmoid fits to subject prediction data';
legend([p1,p2,f1,f2], {txtCond1,txtCond2,txtSigFit1,txtSigFit2},'Location','southeast');
set(gca,'FontName','Helvetica');
hold off;

%% Statistics
%   for each condition, get create summary statistics on the coefficients
%   compare the statistics from the two conditions
%   i.e. 
%       1. First compare slopes
%       2. If the slopes are different, then no need to compare intercepts
%       3. If the slopes are not different (i.e. parallel), then the curves could be 
%           distinguished with distinct intercepts
%

% get rid of null data from excluded subjects
%postParams(postParams(:,1)==0 & postParams(:,2)==0) = [];
%predParams(predParams(:,1)==0 & predParams(:,2)==0) = [];

% Comparing slopes
[slpP, slpH, slpStats] = signrank(postParams(:,2), predParams(:,2));

if slpH == 1
    disp(['Slopes are significantly different with p-value ' num2str(slpP)]);
    slpStats
else
    % Comparing intercepts (useless if slopes are significantly different)
    [intcptP, intcptH, intcptStats] = signrank(postParams(:,1), predParams(:,1));
end

%% violin plots of slopes
figure(2);
hold on;
violinplot([postParams(:,2) predParams(:,2)],'Slope parameter','ViolinColor',[.4 .4 .4],'EdgeColor',[1 1 1],'ShowMean',true);


%% Local Functions
function binNum = repeatBin(signllr)
    % This function assigns the bin number depending on the seqllr value signed by
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

