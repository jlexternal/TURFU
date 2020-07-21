% Testing Confidence models

% Files to import from base workspace when called from another function
expe = evalin('base', 'expe');

% Global variables
global blockFilter;
global globCond;

% For testing purposes %%%%%%%%% DELETE BEFORE ACTUAL USAGE %%%%%%%%%%
globCond = 2;
blockFilter = 'condition';
%%%%%%%%%%%%%%%%%% or your analysis will be shit %%%%%%%%%%%%%%%%%%%%%%%%


% Static variables
condition   = globCond; % 1-postdiction or 2-prediction condition (blockFilter == condition)
                        % 1-low volatility or 2- high volatility (blockFilter == volatility)
n_trials = 72;
haz     = 0.3697; % Hazard rate fitting parameter
sigma   = 0.2452; % Inference noise fitting parameter


% Reference for condition
if strcmpi(blockFilter, 'condition')
    blocks = find([expe.blck.taskid] == condition);	% taskid refers to subject role in exp.
elseif strcmpi(blockFilter, 'volatility')
    blocks = find([expe.blck.condtn] == condition);	% condtn refers to volatility of block
elseif strcmpi(blockFilter, 'all')
    blocks = [3:10];    % when analyzing all non-training blocks
else
    error("Choose appropriate blocking variable: 'condition' (inference direction) or 'volatility'")
end
blocks(find(blocks==1)) = [];   % disregard training blocks
blocks(find(blocks==2)) = [];   % disregard training blocks

% Implementation of the prior function 
phi     = @(lnm1, haz) lnm1 + log((1-haz)./haz + exp(-lnm1)) - log((1-haz)./haz + exp(lnm1));

%% Regular Glaze (2015) hard cut by result from logistic regression

L = zeros(numel(blocks),72);
C = zeros(numel(blocks),72);
hardC = ones(numel(blocks),72);
C2 = zeros(numel(blocks),72);
confid = [];
LLRs = [];

iblock = 1;

drugoctr = 0;
hardctr = 0;
signedctr = 0;
for block = blocks
     % experiment and subject response data for block i
    LLR     = expe.blck(block).seqllr;       % llrs
    LLRs    = horzcat(LLRs, LLR);            
    seqdir  = expe.blck(block).seqdir;       % sequence direction
    rspDir  = expe.rslt(block).resp(2:73);   % subject choice direction: 1-pos, 2-neg
    
    for i = 1:n_trials
        % calculate number of samples in this trial
        n_samples = numel(expe.blck(block).seqtlt{i});
        
        % calculate prior and confidence
        if i == 1
            L(iblock,i) = LLR(1);
            C(iblock,i) = conf(0,n_samples,LLR(i),seqdir(i),sigma);
            C2(iblock,i) = conf2(0,rspDir(i),cutoff);
        else  
            L(iblock,i) = phi(L(iblock,i-1),haz) + LLR(i);
            C(iblock,i) = conf(L(iblock,i-1),n_samples,LLR(i),seqdir(i),sigma);
            C2(iblock,i) = conf2(L(iblock,i-1),rspDir(i),cutoff);
        end
        if abs(L(iblock,i)) > 1.929
            hardC(iblock,i) = hardC(iblock,i)+1;
        end
    
    end
    figure(iblock);
    %scatter([1:72],C(iblock,:),'x');
    scatter([1:72],hardC(iblock,:), '+');
    hold on;
    scatter([1:72],expe.rslt(blocks(iblock)).conf(2:73));
    scatter([1:72],((C2(iblock,:)-1.5).*1.1)+1.5, 'd');
    hold off;
    
drugoctr = drugoctr + numel(find(C(iblock,:)-expe.rslt(blocks(iblock)).conf(2:73)==0))
hardctr = hardctr + numel(find(hardC(iblock,:)-expe.rslt(blocks(iblock)).conf(2:73)==0))
signedctr = signedctr + numel(find(C2(iblock,:)-expe.rslt(blocks(iblock)).conf(2:73)==0))
    
    confid = vertcat(confid,C);
    iblock = iblock+1;
end


%% Confidence on Belief signed by Choice, hard cut by 50% on signed L

L = zeros(numel(blocks),72);

iblock = 1;
for block = blocks
     % experiment and subject response data for block
    LLR     = expe.blck(block).seqllr;       % llrs
    rspdir  = expe.rslt(block).resp(2:73);   % subject choice direction: 1-pos, 2-neg
    
    rspdir  = (rspdir-1.5).*-2;          % convert subj choice: pos -> +1, neg -> -1
    
    for itrial = 1:n_trials
        % calculate posterior (Glaze) signed by choice
        if itrial == 1
            L(iblock,itrial) = LLR(1) .* rspdir(itrial);
        else
            L(iblock,itrial) = phi(L(iblock,itrial-1),haz) + LLR(itrial);
            L(iblock,itrial) = L(iblock,itrial).* rspdir(itrial);
        end
    end
    
    % plotting verification
     
    figure(iblock);
    hold on;
    scatter(L(iblock,:),expe.rslt(blocks(iblock)).conf(2:73));
    hold off;
    
    
    iblock = iblock+1;
end




%% Logistic reg on confidence given signed posterior

L_ccted = [];
C_ccted = [];
LLR_ccted = [];
for i = 1:size(L,1);
    LLR_ccted = horzcat(LLR_ccted, expe.blck(blocks(i)).seqllr);
    L_ccted = horzcat(L_ccted, L(i,:)); % concatenate L values into one array
    C_ccted = horzcat(C_ccted, expe.rslt(blocks(i)).conf(2:73));
end
A = [-6:.1:6];
C_ccted = (C_ccted-1).';
a = glmfit(L_ccted, C_ccted,'binomial','link','logit');
Z = a(1)+(a(2)*A);
Z = 1 ./ (1 + exp(-Z));
cutoff = -a(1)./a(2);

% visual check
clf;
figure(1);
scatter(L_ccted, C_ccted);
hold on;
plot(A,Z);
hold off;

%{
LLR_ccted = abs(LLR_ccted);
a = glmfit(LLR_ccted, C_ccted,'binomial','link','logit');
Z = a(1)+(a(2)*A);
Z = 1 ./ (1 + exp(-Z));
figure(2);
scatter(LLR_ccted, C_ccted);
hold on;
plot(A,Z);
hold off;
%}

%% Creating confidence decision function

% 1. Reflect the logistic function about the y-axis 
% 2. Decision function is the max of the two functions along the support
% 3. Find the minima of the function and find points ±S (symmetric) about the minima
%    that correspond to 50% of choices
% 4. Then, belief that falls between [-S,+S] corresponds to low confidence

function c = conf2(posterior, subjDir, cutoff)
    subjDir = (subjDir-1.5).*-2; % convert subj choice to ±1
    if posterior.*subjDir >= cutoff
        c = 2;
    else
        c = 1;
    end
end
%% Local functions
function c = conf(posterior, n_samples, llr, seqdir,sigm)
    pd  = makedist('Normal', 'mu', posterior, 'sigma', 0); % initialize standard normal distrib.
    
    llr = llr/sqrt(n_samples+sigm.^(-2));  % scale argument of pdf by sqrt(t); t = number of samples*unit time
    
    % deal with high absolute values == favorable belief
    if seqdir == 2
        llr = llr.*-1;
    end
    % calculate cdf (normal) based on entered parameters
    if cdf(pd, llr) > .5 
        c = 2; % high confidence
    else
        c = 1; % low confidence
    end
end