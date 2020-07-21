% Code for Weighted (on the prior/likelihood itself) Bayes LLR Model
%   
%   (NB: The weighting on this model is fundamentally different than weightedBayes2)
%
%   Name:   weightedBayes
%   Type:   Function
% Author:   Jun Seok Lee
%   Date:   February 2019

% ** PREREQUISITES FOR PROPER FUNCTIONALITY OF CODE **
%    Experimental and subject data should be loaded into the BASE workspace

function objFn = weightedBayes(params)

% Files to import from base workspace when called from another function
blck = evalin('base', 'blck');
rslt = evalin('base', 'rslt');
global globCond;    % 'globCond' is a global variable set in a higher-level script
global blockFilter;

% Fitting parameters
w = params;        % weight on prior (belief)

% Initialize static variables
condition   = globCond; % 1-observer/backward or 2-agent/forward condition (blockFilter == condition)
                        % 1-low volatility or 2- high volatility (blockFilter == volatility)
n_trials = 72;
                        
if blockFilter == 'condition'
    blocks = find([blck.taskid] == condition);	% taskid refers to subject role in exp.
elseif blockFilter == 'volatility'
    blocks = find([blck.condtn] == condition);	% condtn refers to volatility of block
else
    err("Choose appropriate blocking variable: 'condition' or 'volatility'")
end

% Objective function to determine best-fitting parameter
% Objective function is a count of non-matching trials to subject data whose
% minimization corresponds to the best fitting model
objFn = 0;

L_rec = [];
% Algorithm:
for i = blocks
    % experiment and subject response data for block i
    LLR     = blck(i).seqllr;       % llrs
    rspDir  = rslt(i).resp(2:73);   % subject choice direction: 1-pos, 2-neg
    % manipulated experiment and subject response data for block i
    signRspDir = ((rspDir-1.5).*-1)*2;  % sign of response direction: 1:pos, -1:neg
    
    objFnInit = 0;
    L = zeros(1,n_trials);
    for j = 1:n_trials
        if j == 1
            L(j) = LLR(j); % first trial belief is the evidence
        else
            L(j) = (w.*L(j-1))+((1-w).*LLR(j)); % Weighted Bayes on prior or likelihood distrib
        end
        
        % Add to objective function if model response does not match human 
        if sign(L(j)) == signRspDir(j)
            objFnInit = objFnInit + 1;
        end
    end
    objFnInit = objFnInit/n_trials;
    objFnInit = log(objFnInit); % test
    objFn = objFn + objFnInit;  % test
    L_rec = vertcat(L_rec, L);  
    
end

assignin('base', 'L_rec', L_rec);
objFn = objFn.*-1;

% Optional visuals for running with BADS
visuals = true; % set to false if you don't want to see the search progression
hold on;
txt = evalin('base', 'txt');
xlim([0 1]);
xlabel('prior weight');
ylabel('negative log likelihood of the model params');
x1 = evalin('base', 'x1');
scatter(w, objFn, '.');
hold on;
text(w, objFn, num2str(txt));
drawnow;
assignin('base', 'txt', txt+1);
hold off;

end