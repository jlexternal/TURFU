% Behavior Simulation for the Future Condition from Parameters found
%   based on the Glaze (2015) model + Inference Noise + Selection Noise

%   Name:   selNoiseGlazeForward_simulator.m
% Author:   Jun Seok Lee
%   Date:   April 2019
%   Type:   Function
%   Note:   For the TURFU experiment

function responses = selNoiseGlazeForward_simulator(haz, sigm, eta, seqllrs, nsamples)

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
            beliefs(itrial)     = seqllrs(itrial) + infNoise(nsamples(itrial)); % L_1
        % for the other trials, update the belief based on the model    
        else
            beliefs(itrial)     = seqllrs(itrial) + psi(beliefs(itrial-1)) + infNoise(nsamples(itrial)); % L_n
        end
        
        %disp(['seqllr: ' num2str(seqllrs(itrial))]); % DEBUG
        %disp(['L_n: ' num2str(beliefs(itrial))]);    % DEBUG
        
        % Forward model (selection noise is introduced here)
        beliefs(itrial) = seqllrs(itrial) + psi(beliefs(itrial)) + selNoise; %L_{n+1} with
        
        %disp(['L_{n+1}: ' num2str(beliefs(itrial))]); % DEBUG
        %pause(1);                                     % DEBUG
        
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
