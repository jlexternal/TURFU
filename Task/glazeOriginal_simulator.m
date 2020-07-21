% Behavior Simulation from Parameters
%   based on the Glaze (2015) model 
%
%   Name:   glazeOriginal_simulator.m
% Author:   Jun Seok Lee
%   Date:   April 2019
%   Type:   Function
%   Note:   Adapted for TURFU experiment
%           Simulates behavior for a single block 

function [responses, beliefs] = glazeOriginal_simulator(haz, seqllrs, futurecond)

% Arguments:
%   haz     := [single/double]  hazard rate parameter
%   seqllrs := [array(single)]  sequence LLRs for the block
    debugging = false; % set to true if debugging

    % Model functions
    psi     = @(lnm1) lnm1 + log((1-haz)./haz + exp(-lnm1)) - log((1-haz)./haz + exp(lnm1));

    n_trials = numel(seqllrs); % number of trials in a given block
    if n_trials ~= 72
       error('Double check the number of values in the seqllr argument!'); 
    end

    % data structures for beliefs and responses
    beliefs   = zeros(1,n_trials);    % belief values (L) from the model simulation
    responses = zeros(1,n_trials+1);  % responses (1/2 for orange/blue) from the model simulation 
    
    responses(1) = 1; % simulating unimportant choice for first answer
    if debugging
        figure(1);
        disp(['Hazard rate: ' num2str(haz)]);
        disp('');
    end

    for itrial = 1:n_trials   
        % first trial belief update
        if itrial == 1
            beliefs(itrial)     = seqllrs(itrial);
            if futurecond
                beliefs(itrial) = psi(beliefs(itrial));
            end
            responses(itrial+1) = belief2decision(beliefs(itrial));  
            continue;
        end
        % for the other trials, update the belief based on the model
        beliefs(itrial)     = psi(beliefs(itrial-1)) + seqllrs(itrial);
        
        if futurecond
            beliefs(itrial) = psi(beliefs(itrial));
        end

        % convert the belief into a decision (1 or 2)
        responses(itrial+1) = belief2decision(beliefs(itrial));    
        
        if debugging 
            disp(num2str(beliefs(itrial)));
            scatter(itrial+1, (responses(itrial+1)-1.5).*-2,'b');
            hold on;
            plot(itrial+1, beliefs(itrial),'Color','r');
            text(itrial+1, beliefs(itrial), num2str(itrial+1));
            xlim([0 73]);
            ylim([-7 7]);
            waitforbuttonpress;
        end
    end
    if debugging
        hold off;
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
