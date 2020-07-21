% Code for Weighted (on the prior/likelihood itself) Bayes LLR Model
%   
%   (NB: The weighting on this model is fundamentally different than weightedBayes2)
%
%   Name:   glazeOriginal
%   Type:   Function
% Author:   Jun Seok Lee
%   Date:   April 2019
%
%   NOTE: This code is modified for the data structures created by the TURFU experiment

% ** PREREQUISITES FOR PROPER FUNCTIONALITY OF CODE **
%    Experimental data should be loaded into the BASE workspace

function objFn = glazeOriginal(params)

% Files to import from base workspace when called from another function
expe        = evalin('base', 'expe');
blockFilter = evalin('base','blockFilter'); % is set to either 'condition' or 'volatility'
globCond    = evalin('base', 'globCond');   % 'globCond' is a global variable set in a higher-level script

% Initialize static variables
condition   = globCond; % 1-postdictive or 2-predictive condition (blockFilter == direction)
                        % 1-low volatility or 2- high volatility (blockFilter == volatility)
n_trials = 72;

% Fitting parameters
haz = params;        % weight on prior (belief)

% Choosing the blocks to analyze
if strcmpi(blockFilter, 'direction')
    blocks = find([expe.blck.taskid] == condition);	% taskid refers to subject role in exp.
elseif strcmpi(blockFilter, 'volatility')
    blocks = find([expe.blck.condtn] == condition);	% condtn refers to volatility of block
else
    error("Choose appropriate blocking variable: 'direction' or 'volatility'")
end
blocks(find(blocks==1)) = [];   % disregard training blocks
blocks(find(blocks==2)) = [];   % disregard training blocks

% Implementation of the phi function 
phi = @(lnm1, haz) lnm1 + log((1-haz)./haz + exp(-lnm1)) - log((1-haz)./haz + exp(lnm1));

% Objective function to determine best-fitting parameter
% Objective function is a count of matching subject trials to the model
objFn = 0;

L_rec = [];
% Algorithm:
for i = blocks
    % experiment and subject response data for block i
    LLR     = expe.blck(i).seqllr;       % llrs
    rspDir  = expe.rslt(i).resp(2:73);   % subject choice direction: 1-pos, 2-neg
    % manipulated experiment and subject response data for block i
    signRspDir = ((rspDir-1.5).*-1)*2;  % sign of response direction: 1:pos, -1:neg
    
    objFnInit = 0;
    L = zeros(1,n_trials);
    for j = 1:n_trials
        if j == 1
            L(j) = LLR(j); % first trial belief is the evidence
        else
            L(j) = phi(L(j-1),haz)+LLR(j); % Original Glaze model (w/o noise)
        end
        
        % Add to objective function if model response matches human response
        if sign(L(j)) == signRspDir(j)
            objFnInit = objFnInit + 1;
        end
    end
    % convert likelihood of data from model into log-likelihood and add for each block
    objFn = objFn + log(objFnInit./72);
    L_rec = vertcat(L_rec, L);  
    
end

assignin('base', 'L_rec', L_rec);
objFn = objFn.*-1; % convert to negative LL 

% Optional visuals for running with BADS
visuals = true; % set to false if you don't want to see the search progression

if visuals == true
    hold on;
    txt = evalin('base', 'txt');
    xlim([0 1]);
    xlabel('hazard rate');
    ylabel('negative objective fn value');
    x1 = evalin('base', 'x1');
    scatter(haz, objFn, '.');
    hold on;
    text(haz, objFn, num2str(txt));
    drawnow;
    assignin('base', 'txt', txt+1);
    hold off;
end

end