function objFn = infNoiseGlaze_PF(params)
% Sequential Monte Carlo Particle Filter on Modified Glaze (2015)
%
%   Name:   infNoiseGlaze_PF
%   Type:   Function
% Author:   Jun Seok Lee
%   Date:   March 2019

% Description:  Generates behavior based on the Glaze (2015) model + Inference Noise
% 
% Necessary data:
% Experiment and subject data
%
% The particle filter, in words:
%   For each trial j 
%       Simulate L_j nOfParticles times for all particles in L_history + replenish(?)
%       If the absolute value of the L_j particle is higher than support limit
%           set its absolute value to 6.5 
%       If the sign of the L_j of particle doesn't match subject choice
%           remove these particles
%           store values in L_history_j
%       If more than 10% of original number of particles are rejected
%           Replenish distribution by factor that normalizes to the original (or lesser)
%            amount to the original number of particles
%       If more than 90% of original number of particles are rejected
%           Consider the subject choice as uninformative
%           Keep original particles before filtering
%
%       ~Loop

% Files to import from base workspace when called from another function
expe = evalin('base', 'expe');

% Global variables
global blockFilter;
global globCond;
global globParticleCount;

% Static Parameters:
n_trials     = 72;              % number of trials per block
n_particles  = globParticleCount;% # of values in the support to run the model on
precision_l  = 0.1;             % grain of the support of L
replenCutoff = 0.9;             % percent*100 threshold of accepted particles before replenishment
uninfoCutoff = 0.1;             % percent*100 threshold of accepted particles triggering uninformative choice
part_thres   = ceil(n_particles .* replenCutoff); % threshold number on accepted particles
getmodelLL  = false;
futurecond = false;

% Reference for condition
if strcmpi(blockFilter, 'direction')
    blocks = find([expe.blck.taskid] == globCond);	% taskid refers to inference direction in exp.
elseif strcmpi(blockFilter, 'volatility')
    blocks = find([expe.blck.condtn] == globCond);	% condtn refers to volatility of block
elseif strcmpi(blockFilter, 'all')
    blocks = [3:10];    % when analyzing all non-training blocks
else
    error("Choose appropriate blocking variable: 'direction' or 'volatility'")
end
blocks(find(blocks==1)) = [];   % disregard training blocks
blocks(find(blocks==2)) = [];   % disregard training blocks

n_blocks = numel(blocks);

hazRate = params(1);            % Modified for BADS
sigma = params(2);              % Modified for BADS

% Optional visuals for running with BADS
if ~getmodelLL
    hold on;
    txt = evalin('base', 'txt');
    xlim([0 1]);
    ylim([0 1]);
    ylabel('sigma');
    xlabel('hazard rate');
    x1 = evalin('base', 'x1');
    if txt ~= 1
        x2 = [hazRate sigma];
        v = x2-x1;
        quiver(x1(1),x1(2),v(1),v(2),0);
    end
    assignin('base', 'x1', [hazRate sigma]);
    scatter(hazRate, sigma, '.');
    hold on;
    text(hazRate, sigma, num2str(txt));
    drawnow;
    assignin('base', 'txt', txt+1);
    hold off;
end


% Objective Function / Output for parameter search
objFn = 0;

% The simulator model:
%   L_i = psi(L_{j-1}, H) + LLR_{n} + eps(m_j)
%     where
%       psi : log-likelihood ratio of prior belief (L_{j-1}) weighted by H
%       LLR : log-likelihood ratio of evidence (defined in loop)
%       err : Gaussian white noise whose stDev = sigma*sqrt(m_j)
%           := random draw from Gaussian(0, sigma*sqrt(m_j))
psi = @(lnm1, haz) lnm1 + log((1-haz)./haz + exp(-lnm1)) - log((1-haz)./haz + exp(lnm1));
err = @(sigm,nsamples) random('Normal', 0, sigm*sqrt(nsamples));

% Data storage and structures (for reliability tracking?):
%   L_acc       = [ {# of accepted particles at each trial} ] 
%                   : array of 72 elements, each element/cell containing variable # of values
%   L_rej       = [ {# of rejected particles at each trial} ]
%                   : array of 72 elements, each element/cell containing variable # of values
%   uninfoCtr   = [ amount of times the subject made an "uninformative choice" ]
%                   : array of 8 elements, each element/cell a counter of "uninformative choice"
L_acc   = {};
L_rej   = {};
uninfoCtr = zeros(1, n_blocks);

% Support and precision of L:
%   The support of L will be [-6.5, 6.5]
%   The grain of the support of L will be determined by lPrecision 
min_L   = -6.5;
max_L   = 6.5;
L_supp  = [min_L:precision_l:max_L];    % binning structure will be [lbound,ubound) 


debugResamplingBool = false;  % resampling error debug message trigger
blckCtr = 1;
% Algorithm:
for i = blocks
    n_part_temp  = n_particles; % block loop starts w/ preset number of particles
    L_acc_btemp  = [];  % temp structure for accepted particles on each block i
    L_rej_btemp  = [];  % temp structure for rejected particles on each block i
    %ll_rej_btemp = 0;   % temp value holder for ll(rejections) on each block
    
    uninformative = false;
    
    %disp(['Simulating block ' num2str(i)]);
    objFn_temp = 0;
    for j = 1:n_trials
        % get subject information
        subjDir = expe.rslt(i).resp(j+1);        % choice of subject; 1-pos, 2-neg
        posDirBool = false;
        
        if subjDir == 1
            posDirBool = true;
        end
        nSmpls  = numel(expe.blck(i).seqang{j}); % # of samples in trial j
        LLR     = expe.blck(i).seqllr(j);        % LLR of evidence in trial j
        
      % Simulation of L's for however many particles is set
        L_acc_ttemp = [];  % temp structure for accepted particles on each trial j
        L_rej_ttemp = [];  % temp structure for rejected particles on each trial j
        L_ttemp     = [];  % temp structure for unfiltered particles on each trial j
        
        % Inference Model (Glaze + noise on llr) **NOT A FORWARD MODEL YET**
        for k = 1:n_part_temp
            if j == 1   % first instance of L has no psi
                L = LLR + err(sigma, nSmpls);
                if futurecond %added for model BIC calc in modelLL_get.m
                    L = psi(L,hazRate);
                end
            else 
                L = LLR + err(sigma, nSmpls) + psi(L_acc_btemp{1,j-1}(k), hazRate);
                if futurecond %added for model BIC calc in modelLL_get.m
                    L = psi(L,hazRate);
                end
            end
            L_ttemp = horzcat(L_ttemp, L);
            
            % Filter particles that correspond to the subject choice 
            if (posDirBool == true & L > 0) | (posDirBool == false & L < 0)
                L_acc_ttemp = horzcat(L_acc_ttemp, L);
            else
                L_rej_ttemp = horzcat(L_rej_ttemp, L);
            end
            
        end
        
        % Trigger for when too many particles were rejected due to uninformative choice
        if numel(L_acc_ttemp) < uninfoCutoff.*n_particles
            uninformative = true;
%            disp(['  ' num2str(numel(L_acc_ttemp)) ' particles were accepted.']);
            uninfoCtr(blckCtr) = uninfoCtr(blckCtr) + 1;    % keep track of how many times this was triggered for a subject
%            disp(['  Uninformative choice made at block ' num2str(i) ', trial ' num2str(j)]);
            L_acc_ttemp = horzcat(L_acc_ttemp, L_rej_ttemp);    % keep all rejected particles
%            disp(['  ' num2str(numel(L_acc_ttemp)) ' are being passed through.']);
            L_rej_ttemp = [];                                   % reject no particles
        end
        
%         %%%%%%%% DEBUG code: A check to see if 0 particles were accepted %%%%%%%%
%         if numel(L_acc_ttemp) == 0                                       %%%%%%%%
%             disp(['block ' num2str(i) ' trial ' num2str(j)]);            %%%%%%%%
%             disp(['Number of particles at trial ' num2str(j-1) ' was ' num2str(numel(L_acc_btemp{1, j-1}))]);
%             error(['Error: No particles were accepted. Possibly due to low particle count' ...
%                   ' in the previous trial.']);                           %%%%%%%%
%         end                                                              %%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%% Delete when fixed %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Objective Function Update
        %   Adds the exponents of the fraction of total particles accepted to be
        %   converted to log later for LogSumExp
        if ~uninformative
            n_part_acc = numel(L_acc_ttemp);
            if ~getmodelLL
                objFn_temp = objFn_temp + exp(n_part_acc./n_part_temp); % exp used for LogExpSum
            else
                objFn_temp = objFn_temp + log(n_part_acc./n_part_temp); % standard log-likelihood
            end
        end
            
        % Replenishing condition for particles
        % **replenished particles hold random values within bin limits [lbnd,ubnd) **
        nFiltered = numel(L_acc_ttemp);
        if nFiltered < part_thres               % if number of accepted particles less than threshold
            replenFactor = n_particles./nFiltered;   % calculate replenishing factor
            
%             disp(['block ' num2str(i) ' trial ' num2str(j)]);                   %DEBUG TEST CODE
%             disp('Accepted number of particles is less than threshold');        %DEBUG TEST CODE
%             disp(['current number of particles is ' num2str(nFiltered)]);       %DEBUG TEST CODE
%             disp(['replenishing factor: ' num2str(rFactor)]);                   %DEBUG TEST CODE
            
            L_acc_ttemp = sort(L_acc_ttemp);      % sort the particles for running replenishment
            
            % running replenishment after the particles are ordered
            replen_temp = [];
            
            % Bin range considerations
            if posDirBool == false  % if negative choice, support less than 0
                suppStartPt = 1;
                suppEndPt   = numel(L_supp) - ceil(numel(L_supp)./2);
            else                    % if positive choice, support greater than 0
                suppStartPt = ceil(numel(L_supp)./2);
                suppEndPt   = numel(L_supp);
            end
            
            % Given that the subject made an uninformative choice, override above,
            %  and examine full support
            if uninformative 
                suppStartPt = 1;
                suppEndPt   = numel(L_supp);
                uninformative = false;
            end
            
            partCtr     = 0;    % counts number of particles in a given bin
            partPointer = 1;    % points to position of particle in filtered set
            for k = suppStartPt:suppEndPt-1 % loop through the appropriate bins
                lBound = L_supp(k);   % lower inclusive bound of bin [
                uBound = L_supp(k+1); % upper exclusive bound of bin )
                replenish = false;
                while ~replenish 
                    if partPointer == nFiltered
                        break
                    end
                    if L_acc_ttemp(partPointer) >= lBound && L_acc_ttemp(partPointer) < uBound
                        partCtr     = partCtr + 1;
                        partPointer = partPointer + 1;
                    else
                        replenish = true;
                        break
                    end
                end
                if replenish
                    discrep = (replenFactor.*partCtr) - partCtr; % how much replenishment to normalize
                    discrep = round(discrep); % round this number to integer
                    if discrep > 0
                        replen_temp = horzcat(replen_temp, unifrnd(lBound,uBound,1,discrep));
                    end
                end
                partCtr = 0;
            end
            
           %%%%%%%% DEBUG code: A check to see if 0 particles were accepted %%%%%%%%%%%%%%%
             if debugResamplingBool                                                  %%%%%%
                 disp(['number before replenishment ' num2str(numel(L_acc_ttemp))]); %%%%%%
             end                                                                     %%%%%%

            L_acc_ttemp = horzcat(L_acc_ttemp, replen_temp); % DO NOT DELETE w/ DEBUG CODE
            
             if debugResamplingBool                                                  %%%%%%
                 disp(['number after replenishment ' num2str(numel(L_acc_ttemp))]);  %%%%%%
             end                                                                     %%%%%%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%% Delete when fixed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        L_ttemp      = L_acc_ttemp;
        L_acc_btemp  = horzcat(L_acc_btemp, {L_acc_ttemp});
        n_part_temp  = numel(L_acc_ttemp);
        L_rej_btemp  = horzcat(L_rej_btemp, {L_rej_ttemp});
        
         %%%%%%%% DEBUG code: A check to see if small amount of particles accepted %%%%%%%%%%%%%%%
%         if n_part_temp < 0.1.*n_particles
%             disp(['number of passed particles are ' num2str(n_part_temp) ' at block ' num2str(i) ' trial ' num2str(j)]);
%             debugResamplingBool = true;
%         else
%             debugResamplingBool = false;
%         end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%% Delete when fixed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    if ~getmodelLL
        objFn = objFn + log(objFn_temp);    % convert to LogSumExp
    else
        objFn = objFn + objFn_temp;         % keep as log-likelihood
    end
    %disp(['objFn after block ' num2str(objFn)]);
    
    L_acc = vertcat(L_acc, L_acc_btemp);
    L_rej = vertcat(L_rej, L_rej_btemp);
    blckCtr = blckCtr + 1;
end

percUninf = sum(uninfoCtr)/(n_blocks.*n_trials).*100;
%disp([num2str(sum(uninfoCtr)) ' uninformative trials were found.']);
%disp([num2str(percUninf) ' percent of choice across all blocks uninformative.']);

% Files to export to base workspace when called from another function
assignin('base', 'L_acc', L_acc);
assignin('base', 'L_rej', L_rej);
assignin('base','percUninf',percUninf);

% Objective function output
if objFn == -inf | objFn == inf
    objFn = 0;
else
    objFn = objFn.*-1; % convert to negative for minimum search on BADS
end

%disp(num2str(objFn));

end