% Code for Weighted (on the prior/likelihood itself) Bayes LLR Model
%   
%   (NB: The weighting on this model is fundamentally different than weightedBayes2)
%
%   Name:   glazeOriginal_v_modelRec
%   Type:   Function
% Author:   Jun Seok Lee
%   Date:   April 2019
%
%   NOTE: This code is modified for the data structures created by the TURFU experiment

function objFn = glazeOriginal_v_modelRec(params)

% Files to import from base workspace when called from another function
simBehavior = evalin('base', 'simBehavior');
hIndex      = evalin('base', 'hIndex');
expeIndex   = evalin('base', 'expeIndex');

n_trials = 72;

% Fitting parameters
haz = params;        % weight on prior (belief)

% Implementation of the phi function 
psi = @(lnm1, haz) lnm1 + log((1-haz)./haz + exp(-lnm1)) - log((1-haz)./haz + exp(lnm1));

% Objective function to determine best-fitting parameter
% Objective function is a count of non-matching trials to subject data whose
% minimization corresponds to the best fitting model
objFn = 0;

L_rec = [];
% Algorithm:
% experiment and subject response data 
LLR     = simBehavior(hIndex).expe(expeIndex).seqllr; % llrs
rspDir  = simBehavior(hIndex).rslt(expeIndex).resp;   % subject choice direction: 1-pos, 2-neg

% manipulated experiment and subject response data for block i
signRspDir = ((rspDir-1.5).*-1)*2;  % sign of response direction: 1:pos, -1:neg

objFnInit = 0;
L = zeros(1,n_trials);
for j = 1:n_trials
    if j == 1
        L(j) = LLR(j); % first trial belief is the evidence
    else
        L(j) = psi(L(j-1),haz)+LLR(j); % Original Glaze model (w/o noise)
    end

    % Add to objective function if model response does not match human 
    if sign(L(j)) == signRspDir(j+1)
        objFnInit = objFnInit + 1;
    end
end

%objFnInit = objFnInit/n_trials;
%objFnInit = log(objFnInit); % test
objFn = objFnInit;  % test
L_rec = vertcat(L_rec, L);  
   

assignin('base', 'L_rec', L_rec);
objFn = objFn.*-1;

% Optional visuals for running with BADS
visuals = true; % set to false if you don't want to see the search progression

if visuals == true
    hold on;
    txt = evalin('base', 'txt');
    xlim([0 1]);
    xlabel('hazard rate');
    ylabel('negative log likelihood of the model params');
    x1 = evalin('base', 'x1');
    scatter(haz, objFn, '.');
    hold on;
    text(haz, objFn, num2str(txt));
    drawnow;
    assignin('base', 'txt', txt+1);
    hold off;
end
end