% Postdictive inference parameter search for a single subject

%   Name:   postdictive.m
%   Type:   Script
% Author:   Jun Seok Lee
%   Date:   March 2019

% 1. Load up subject
% 2. Set inference direction
% 3. Run DETERMINISTIC Glaze (2015) model to get hazard rate (H)
% 4. Logistic regression on confidence with posteriors (L) given by the model as predictors
%   4a. Predictor as absolute value of L
%   4b. Predictor as L signed with subject choice 
% 5. Compare 4a and 4b

%% 1. Load experiment and subject data:
subject = 6; % choose from 1 through #
global globSubj;
global globCond;
global blockFilter;
globSubj = subject;
globCond = 1;
n_trials = 72;

%% 2. Set inference direction
blockFilter = 'direction';
modelRec = false;

%% 3. Get hazard rate parameter:
%   This runs glazeOriginal_paramSearchBADS to find the optimalhazard rate parameter 
H = 0;
for i = 1:5
    run glazeOriginal_paramSearchBADS;
    H = H+optiParams;
end
H = H./5;

%% 4. Run logistic regression 
%   4a. with abs(L)                 as predictor variable
%   4b. with L signed with choice   as predictor variable

% Prior function of Glaze et al. (2015)
phi     = @(lnm1, haz) lnm1 + log((1-haz)./haz + exp(-lnm1)) - log((1-haz)./haz + exp(lnm1));

% Choosing the blocks corresponding to analysis
if strcmpi(blockFilter, 'direction')
    blocks = find([expe.blck.taskid] == globCond);	% taskid refers to subject role in exp.
elseif strcmpi(blockFilter, 'volatility')
    blocks = find([expe.blck.condtn] == globCond);	% condtn refers to volatility of block
else
    error("Choose appropriate blocking variable: 'direction' or 'volatility'")
end
blocks(find(blocks==1)) = [];   % disregard training blocks
blocks(find(blocks==2)) = [];   % disregard training blocks

L_4a = zeros(numel(blocks),72);
L_4a_ccted = [];
L_4b = zeros(numel(blocks),72);
L_4b_ccted = [];
C_ccted = [];

blckCtr = 1;
for block = blocks
     % experiment and subject response data for block
    LLR     = expe.blck(block).seqllr;      % llrs
    rspdir  = expe.rslt(block).resp(2:73);  % subject choice direction: 1-pos, 2-neg 
    rspdir  = (rspdir-1.5).*-2; % convert subj choice: pos -> +1, neg -> -1
    
    for itrial = 1:n_trials
        % calculate posterior (Glaze)
        if itrial == 1
            L_4a(blckCtr,itrial) = abs(LLR(1));
            L_4b(blckCtr,itrial) = LLR(1) .* rspdir(itrial);
        else
            L_4a(blckCtr,itrial) = abs(phi(L_4a(blckCtr,itrial-1),H) + LLR(itrial));
            L_4b(blckCtr,itrial) = phi(L_4b(blckCtr,itrial-1),H) + LLR(itrial);
            L_4b(blckCtr,itrial) = L_4b(blckCtr,itrial).* rspdir(itrial);
        end
    end
    C_ccted = horzcat(C_ccted, expe.rslt(block).conf(2:73));
    L_4a_ccted = horzcat(L_4a_ccted, L_4a(blckCtr,:));
    L_4b_ccted = horzcat(L_4b_ccted, L_4b(blckCtr,:));
    blckCtr = blckCtr+1;
end

xAxis = [-7:.1:7];
C_ccted = (C_ccted-1).';
a = glmfit(L_4a_ccted, C_ccted,'binomial','link','logit');
Z_4a = a(1)+(a(2)*xAxis);
Z_4a = 1 ./ (1 + exp(-Z_4a));
cutoff_4a = -a(1)./a(2);

b = glmfit(L_4b_ccted, C_ccted,'binomial','link','logit');
Z_4b = b(1)+(b(2)*xAxis);
Z_4b = 1 ./ (1 + exp(-Z_4b));
cutoff_4b = -b(1)./b(2);

% visual check
clf;
figure(1);
hold on;
plot(xAxis,Z_4a);
plot(xAxis,Z_4b);
hold off;

%% 5. Confidence model comparison

conf_4a = zeros(numel(C_ccted),1);
conf_4b = zeros(numel(C_ccted),1);

for i = 1:numel(C_ccted)
    % below cutoff scenario
    if L_4a_ccted(i)>cutoff_4a
        conf_4a(i) = 1;
    end
    if L_4b_ccted(i)>cutoff_4b
        conf_4b(i) = 1;
    end
end

disp(['Cutoff based on abs. value of L accounted for ' num2str(numel(find(conf_4a-C_ccted ==0))) ' out of ' num2str(288) ' trials.']);
disp(['Cutoff based on choice signed value of L accounted for ' num2str(numel(find(conf_4b-C_ccted ==0))) ' out of ' num2str(288) ' trials.']);
 