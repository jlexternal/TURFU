function objFn = infNoiseGlaze_PFpar(params)
% Sequential Monte Carlo Particle Filter on Modified Glaze (2015)
% (Code modified from the original for Parallel Processing)
% Author:   Jun Seok Lee
% Date:     February 2019
% 
% Necessary data:
%   Run these commands in the folder containing subject folders 
%       ->  filename = sprintf('S##/expe.mat');
%       ->  load(filename);
%   where ## corresponds to the subject number.

% Files to import from base workspace when called from another function
blck = evalin('base', 'blck');
rslt = evalin('base', 'rslt');

% Static Parameters:
n_blocks     = numel(blck);     % number of blocks per subject
n_trials     = 72;              % number of trials per block
n_particles  = 2000;            % # of values in the support to run the model on
precision_l  = 0.1;             % grain of the support of L
replenCutoff = 0.9;             % percent*100 threshold of accepted particles before replenishment
uninfoCutoff = 0.1;             % percent*100 threshold of accepted particles triggering uninformative choice
part_thres   = ceil(n_particles .* replenCutoff); %threshold number on accepted particles
%hazRate      = 0.3;             % to-be-fitted parameter, H in Glaze (2015) model
%sigma        = 0.4;             % to-be-fitted parameter, subjective noise on llr of evidence

% Modified for BADS
sigma = params(1);
hazRate = params(2);

% Objective Function / Output for parameter search
ll_rej = 0;

% Optional visuals for running with BADS
txt = evalin('base', 'txt');
xlim([0 1]);
ylim([0 1]);
xlabel('sigma');
ylabel('hazard rate');
scatter(sigma, hazRate, '.');
hold on;
text(sigma, hazRate, num2str(txt));
drawnow;
assignin('base', 'txt', txt+1);

% The simulator model:
%   L_i = Phi(L_{j-1}, H) + LLR_{n} + eps(m_j)
%     where
%       phi : log-likelihood ratio of prior belief (L_{j-1}) weighted by H
%       LLR : log-likelihood ratio of evidence (defined in loop)
%       err : Gaussian white noise whose stDev = sigma*sqrt(m_j)
%           := random draw from Gaussian(0, sigma*sqrt(m_j))
phi = @(lnm1, haz) lnm1 + log((1-haz)./haz + exp(-lnm1)) - log((1-haz)./haz + exp(lnm1));
err = @(sigm,nSamples) random('Normal', 0, sigm*sqrt(nSamples));

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

% Support and precision:
%   The support of L will be [-6.5, 6.5]
%   The grain of the support of L will be determined by lPrecision 
min_L   = -6.5;
max_L   = 6.5;
L_supp  = [min_L:precision_l:max_L];    % binning structure will be [lbound,ubound) 

debugResamplingBool = false;  % resampling error debug message trigger

% Algorithm:
for i = 1:n_blocks
    n_part_temp  = n_particles; % block loop starts w/ preset number of particles
    L_acc_btemp  = [];  % temp structure for accepted particles on each block i
    L_rej_btemp  = [];  % temp structure for rejected particles on each block i
    
    uninformative = false;
    
    disp(['Simulating block ' num2str(i)]);
    
    for j = 1:n_trials
        subjDir = rslt(i).resp(j+1);        % choice of subject; 1-pos, 2-neg
        posDirBool = false;
        
        if subjDir == 1
            posDirBool = true;
        end
        nSmpls  = numel(blck(i).seqang{j}); % # of samples in trial j
        LLR     = blck(i).seqllr(j);        % LLR of evidence in trial j
        
      % Simulation of L's for however many particles is set
        par_L_acc_ttemp = zeros(1,n_part_temp);  % temp structure for accepted particles on each trial j
        par_L_rej_ttemp = zeros(1,n_part_temp);  % temp structure for rejected particles on each trial j
        par_L_ttemp = zeros(1,n_part_temp);  % temp structure for unfiltered particles on each trial j
        
        % Inference Model (Glaze + noise on llr) ***NOT A FORWARD MODEL YET***
        %   Parallel Processing Variable I/O
        %       Inputs: [N particle values] from trial j-1              -> par_L_ttemp
        %      Outputs: [N updated particle values] from trial j        
        %               [M accepted particle values] from trial j
        %               [Q rejected particle values] from trial j
        
        par_L = zeros(1,n_part_temp);   % for each trial, create zero array w/ number of unfiltered particles (par. proc.)
%        if j~=1
%            disp(L_acc_btemp{1,j-1});
%        end
        parfor k = 1:n_part_temp
            if j == 1   % first instance of L has no phi
                par_L(k) = LLR + err(sigma, nSmpls);
            else 
                par_L(k) = LLR + err(sigma, nSmpls) + phi(L_acc_btemp{1,j-1}(k),hazRate);
            end
            par_L_ttemp(k) = par_L(k);
            
            % Filter particles that correspond to the subject choice 
            if (posDirBool == true & par_L(k) > 0) | (posDirBool == false & par_L(k) < 0)
                par_L_acc_ttemp(k) = par_L(k);
            else
                par_L_rej_ttemp(k) = par_L(k);
            end
        end
        
        par_L_acc_ttemp = par_L_acc_ttemp(par_L_acc_ttemp~=0);
        par_L_rej_ttemp = par_L_rej_ttemp(par_L_rej_ttemp~=0);
        
        
        % Trigger for when too many particles were rejected due to uninformative choice
        if numel(par_L_acc_ttemp) < uninfoCutoff.*n_particles
            uninformative = true;
            uninfoCtr(i) = uninfoCtr(i) + 1;    % keep track of how many times this was triggered for a subject
            par_L_acc_ttemp = horzcat(par_L_acc_ttemp, par_L_rej_ttemp);    % keep all rejected particles
            par_L_rej_ttemp = [];                                   % reject no particles
%            disp(['  ' num2str(numel(par_L_acc_ttemp)) ' particles were accepted.']);
%            disp(['  Uninformative choice made at block ' num2str(i) ', trial ' num2str(j)]);
%            disp(['  ' num2str(numel(par_L_acc_ttemp)) ' are being passed through.']);
        end
        
        %%%%%%%% DEBUG code: A check to see if 0 particles were accepted %%%%%%%%
%         if numel(par_L_acc_ttemp) == 0                                       %%%%%%%%
%             disp(['block ' num2str(i) ' trial ' num2str(j)]);            %%%%%%%%
%             disp(['Number of particles at trial ' num2str(j-1) ' was ' num2str(numel(L_acc_btemp{1, j-1}))]);
%             error(['Error: No particles were accepted. Possibly due to low particle count' ...
%                   ' in the previous trial.']);                           %%%%%%%%
%         end                                                              %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% Delete when fixed %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Log Likelihood of ~data|parameter / Objective Function Update
        n_part_rej = numel(par_L_rej_ttemp);
        if n_part_rej ~= 0  % in case there are no rejected particles, to avoid log(0)
            ll_rej = ll_rej + log(n_part_rej./n_part_temp);
        end
        
        % Replenishing condition for particles
        % **replenished particles hold random values within bin limits [lbnd,ubnd) **
        nFiltered = numel(par_L_acc_ttemp);
        if nFiltered < part_thres               % if number of accepted particles less than threshold
            rFactor = n_particles./nFiltered;   % calculate replenishing factor
            
%             disp(['block ' num2str(i) ' trial ' num2str(j)]);                   %DEBUG TEST CODE
%             disp('Accepted number of particles is less than threshold');        %DEBUG TEST CODE
%             disp(['current number of particles is ' num2str(nFiltered)]);       %DEBUG TEST CODE
%             disp(['replenishing factor: ' num2str(rFactor)]);                   %DEBUG TEST CODE
            
            par_L_acc_ttemp = sort(par_L_acc_ttemp);      % sort the particles for running replenishment
            
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
                    if par_L_acc_ttemp(partPointer) >= lBound && par_L_acc_ttemp(partPointer) < uBound
                        partCtr     = partCtr + 1;
                        partPointer = partPointer + 1;
                    else
                        replenish = true;
                        break
                    end
                end
                if replenish
                    discrep = (rFactor.*partCtr) - partCtr; % how much replenishment to normalize
                    discrep = round(discrep); % round this number to integer
                    if discrep > 0
                        replen_temp = horzcat(replen_temp, unifrnd(lBound,uBound,1,discrep));
                    end
                end
                partCtr = 0;
            end
            
%           %%%%%%%% DEBUG code: A check to see if 0 particles were accepted %%%%%%%%%%%%%%%
%             if debugResamplingBool                                                  %%%%%%
%                 disp(['number before replenishment ' num2str(numel(par_L_acc_ttemp))]); %%
%             end                                                                     %%%%%%

            par_L_acc_ttemp = horzcat(par_L_acc_ttemp, replen_temp); % DO NOT DELETE w/ DEBUG CODE      
            
%             if debugResamplingBool                                                  %%%%%%
%                 disp(['number after replenishment ' num2str(numel(par_L_acc_ttemp))]);  %%
%             end                                                                     %%%%%%
%           %%%%%%%%%%%%%%%%%%%%%%%%%%%% Delete when fixed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        L_ttemp      = par_L_acc_ttemp;
        L_acc_btemp  = horzcat(L_acc_btemp, {par_L_acc_ttemp});
        n_part_temp  = numel(par_L_acc_ttemp);
        L_rej_btemp  = horzcat(L_rej_btemp, {par_L_rej_ttemp});
        
%         %%%%%%%% DEBUG code: A check to see if small amount of particles accepted %%%%%%%%%%%%%%%
%         if n_part_temp < 0.1.*n_particles    % DEBUG code
%             disp(['number of passed particles are ' num2str(n_part_temp) ' at block ' num2str(i) ' trial ' num2str(j)]);
%             debugResamplingBool = true;
%         else
%             debugResamplingBool = false;
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%% Delete when fixed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    L_acc = vertcat(L_acc, L_acc_btemp);
    L_rej = vertcat(L_rej, L_rej_btemp);
end

percUninf = sum(uninfoCtr)/(n_blocks.*n_trials).*100;
disp(['For this subject, ' num2str(sum(uninfoCtr)) ' uninformative trials were found.']);
disp([num2str(percUninf) ' percent of choice across all blocks uninformative.']);
disp(['All blocks (' num2str(n_blocks) ') simulated.'])


objFn = ll_rej;

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
end