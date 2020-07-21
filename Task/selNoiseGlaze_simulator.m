% Behavior Simulation for the Future Condition from Parameters found
%   based on the Glaze (2015) model + Inference Noise + Selection Noise in the Past Condition
%
%   Name:   selNoiseGlaze_simulator.m
% Author:   Jun Seok Lee
%   Date:   April 2019
%   Type:   Function
%   Note:   For the TURFU experiment

function [responses, beliefs] = selNoiseGlaze_simulator(haz, sigm, eta, seqllrs, nsamples, futurecond)

% Arguments:
%   haz     := [single/double]  hazard rate parameter
%   sigm    := [single/double]  inference noise parameter
%   eta     := [single/double]  selection noise parameter
%   seqllrs := [array(single)]  sequence LLRs for the block
%   nsamples:= [array{integer)] number of samples on a given trial

    % Model functions
    psi      = @(lnm1) lnm1 + log((1-haz)./haz + exp(-lnm1)) - log((1-haz)./haz + exp(lnm1));
    infNoise = @(nsamples) random('Normal', 0, sigm*sqrt(nsamples));
    selNoise = random('Normal', 0, eta);

    n_trials = numel(seqllrs); % number of trials in a given block
    if n_trials ~= 72
       error('Double check the number of values in the seqllr argument!'); 
    end

    % data structures for beliefs and responses
    beliefs   = zeros(1,n_trials);  % belief values (L) from the model simulation
    responses = zeros(1,n_trials+1);  % responses (1/2 for orange/blue) from the model simulation 
    
    responses(1) = 1; % simulating unimportant choice for first answer

    for itrial = 1:n_trials
        % first trial belief update
        if itrial == 1
            beliefs(itrial) = seqllrs(itrial); % L_1
        % for the other trials, update the belief based on the model    
        else
            beliefs(itrial) = seqllrs(itrial) + psi(beliefs(itrial-1)); % L_n
        end
        % add inference and selection noise 
        beliefs(itrial) = beliefs(itrial) + infNoise(nsamples(itrial));
        if futurecond
            beliefs(itrial) = psi(beliefs(itrial));
        end
        beliefs(itrial) = beliefs(itrial) + selNoise; 
        
        % convert the belief into a decision (1 or 2)
        responses(itrial+1) = belief2decision(beliefs(itrial));    
    end
    if numel(responses(responses==0)) ~= 0
        error('Null response(s) found. Check your code.');
    end
    
end

% This function converts a belief value into a binary decision
function decision = belief2decision(belief)
    finished = false;
    while ~finished
        if belief > 0 
            decision = 1;
            finished = true;
        elseif belief < 0 
            decision = 2;
            finished = true;
        else
            belief = rand()-.5;
            continue;
        end
    end
end
