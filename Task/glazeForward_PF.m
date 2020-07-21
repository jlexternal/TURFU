% Forward Inference Model based on Glaze (2015)
% (Particle Filter Simulator)
%
%   Name:  	glazeForward
%   Type: 	Function
% Author:   Jun Seok Lee
%   Date:  	March 2019

%  Model:   L_i = Phi(L_{j-1}, H) + LLR_{n} + eps(m_j)

% ** PREREQUISITES FOR PROPER FUNCTIONALITY OF CODE **
%       Experimental and subject data should be loaded into the BASE workspace

function objFn = glazeForward_PF(params)
% params should be provided by the script running BADS (or something else)

blck = evalin('base', 'blck');
rslt = evalin('base', 'rslt');
subject = evalin('base', 'subject');

% Instantiate global variables declared in higher-level scripts
global blockFilter; % 'blockFilter' is the name of the condition in consideration
global globCond;    % 'globCond' is the numerical value for the condition in question 

% Initialize static variables
condition   = globCond; % 1-postdiction or 2-prediction condition (blockFilter == 'condition')
                        % 1-low or 2- high volatility condition   (blockFilter == 'volatility')
n_trials    = 72;
n_particles = 1000;
sigma       = params(1);        % Modified for BADS
hazRate     = params(2);        % Modified for BADS
uninfoCutoff = 0.1;             % percent*100 threshold of accepted particles triggering uninformative choice

if blockFilter == 'condition'
    blocks = find([blck.taskid] == condition);	% taskid refers to subject role in exp.
elseif blockFilter == 'volatility'
    blocks = find([blck.condtn] == condition);	% condtn refers to volatility of block
else
    error("Choose appropriate blocking variable: 'condition' or 'volatility'");
end

% Anonymous functions:
%       phi : log-likelihood ratio of prior belief (L_{j-1}) weighted by H
%       err : Gaussian white noise whose stDev = sigma*sqrt(m_j) % random draw from Gaussian(0, sigma*sqrt(m_j))
phi = @(lnm1, haz) lnm1 + log((1-haz)./haz + exp(-lnm1)) - log((1-haz)./haz + exp(lnm1));
err = @(sigm, nSamples) random('Normal', 0, sigm*sqrt(nSamples));

% Particle Filter Settings:
%  Replenishment:
replenCutoff    = 0.9; % percent*100 threshold of accepted particles before replenishment
part_thres      = ceil(n_particles .* replenCutoff); % threshold number on accepted particles           
%   The support of L is [-6.5, 6.5]
precision_l = 0.1; % grain of the support of L
min_L       = -6.5;
max_L       =  6.5;
L_supp      = [min_L:precision_l:max_L]; % binning structure will be [lbound,ubound) 

% Sanity check structures
L_acc   = {};
L_rej   = {};
uninfoCtr = zeros(1, numel(blocks));    % tracks number of uninformative choices made by subject

% Optional visuals for running with BADS
hold on;
txt = evalin('base', 'txt');
xlim([0 1]);
ylim([0 1]);
xlabel('sigma');
ylabel('hazard rate');
x1 = evalin('base', 'x1');
if txt ~= 1
    x2 = [sigma hazRate];
    v = x2-x1;
    quiver(x1(1),x1(2),v(1),v(2),0);
end
assignin('base', 'x1', [sigma hazRate]);
scatter(sigma, hazRate, '.');
hold on;
text(sigma, hazRate, num2str(txt));
drawnow;
assignin('base', 'txt', txt+1);
hold off;

% Initiate temporary objective function
ll_rej = 0;

% Algorithm:
% At time step n, the subject makes a postdiction via the inference model à la Glaze 
% about what just happened. Then, it passes the output of that postdiction another
% time into the inference model and uses that output as its prediction.

for i = blocks
    i_index = 0;
    i_index = i_index + 1;
    n_part_temp  = n_particles; % block loop starts w/ preset number of particles
    L_acc_btemp  = [];  % temp structure for accepted particles on each block i
    L_rej_btemp  = [];  % temp structure for rejected particles on each block i
    
    uninformative = false; % presume that subjects have made informative choices
    LLR = blck(i).seqllr;
    
    for j = 1:n_trials
        subjDir = sign(rslt(i).resp(j+1));	% choice of subject; 1-pos, 2-neg
        n_samples = numel(blck(i).seqang{j});
        
        % Simulation of L's for however many particles is set
        L_acc_ttemp = [];  % temp structure for accepted particles on each trial j
        L_rej_ttemp = [];  % temp structure for rejected particles on each trial j
        
        % Generate particles for model simulation
        for k = 1:n_part_temp
            % Simulate the trial for a particle
            if j == 1
                L = LLR(j) + err(sigma, n_samples); % postdiction
                L = phi(L, hazRate) + LLR(j); % prediction
            else
                L = phi(L, hazRate) + LLR(j) + err(sigma, n_samples); % postdiction
                L = phi(L, hazRate) + LLR(j); % prediction
            end
            % Filter accepted and rejected particles based on subject choice
            if sign(L) == subjDir
                L_acc_ttemp = horzcat(L_acc_ttemp, L);  %store in accepted struct
            else
                L_rej_ttemp = horzcat(L_rej_ttemp, L);  %store in rejected struct
            end
        end
        
        % Cneck that a sufficient number of particles were accepted
        if numel(L_acc_ttemp) < uninfoCutoff.*n_particles
            uninformative = true;
            uninfoCtr(i_index) = uninfoCtr(i_index) + 1;    % keep track of how many times this was triggered for a subject
            L_acc_ttemp = horzcat(L_acc_ttemp, L_rej_ttemp);    % keep all rejected particles
            L_rej_ttemp = [];                                   % reject no particles
        end
        
        % Objective Function Update ( log likelihood of p(data|parameter) ) 
        n_part_rej = numel(L_rej_ttemp);
        if n_part_rej ~= 0 & ~uninformative % IN CASE there are no accepted particles, to avoid log(0)
            ll_rej = ll_rej + log(n_part_rej./n_part_temp);
        else % make a file w/ the instance info, but don't break the loop
            fileID = fopen('exp.txt','w');
            fprintf(fileID,'Log 0 encountered for subject %d on block %d, trial %d.',subject,i, j);
            fclose(fileID);
        end
        
        % Replenishment condition
        nFiltered = numel(L_acc_ttemp);
        if nFiltered < part_thres               %"If number of accepted particles less than threshold..."
            rFactor = n_particles./nFiltered;   % calculate replenishing factor
            L_acc_ttemp = sort(L_acc_ttemp);    % sort the particles for running replenishment
            
            % Running replenishment after the particles are ordered
            replen_temp = [];
            % Replenish all over, only on left, or only on right
            if uninformative % examine full support
                suppStartPt = 1;
                suppEndPt   = numel(L_supp);
                uninformative = false;
            else
                if subjDir == -1  % if negative choice, support less than 0
                    suppStartPt = 1;
                    suppEndPt   = numel(L_supp) - ceil(numel(L_supp)./2);
                else                    % if positive choice, support greater than 0
                    suppStartPt = ceil(numel(L_supp)./2);
                    suppEndPt   = numel(L_supp);
                end
            end
            
            partCtr     = 0;    % counts number of particles in a given bin
            partPointer = 1;    % points to position of particle in filtered set
            % The replenishing process occurs BIN by BIN
            for k = suppStartPt:suppEndPt-1 % loop through the appropriate bins
                lBound = L_supp(k);   % lower inclusive bound of bin [
                uBound = L_supp(k+1); % upper exclusive bound of bin )
                replenish = false;    % Boolean "switch" to begin replenishing process
                while ~replenish 
                    if partPointer == nFiltered % "If all particles are accounted for..."
                        break
                    end
                    if L_acc_ttemp(partPointer) >= lBound && L_acc_ttemp(partPointer) < uBound
                        %"If the particle I'm looking at is within my considered bin..."
                        partCtr     = partCtr + 1;
                        partPointer = partPointer + 1;
                    else %"... otherwise, replenish the bin."
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
            
            L_acc_ttemp = horzcat(L_acc_ttemp, replen_temp); 
        end
        
        L_acc_btemp  = horzcat(L_acc_btemp, {L_acc_ttemp});
        n_part_temp  = numel(L_acc_ttemp);
        L_rej_btemp  = horzcat(L_rej_btemp, {L_rej_ttemp});
    end
    L_acc = vertcat(L_acc, L_acc_btemp);
    L_rej = vertcat(L_rej, L_rej_btemp);
    ll_rej = ll_rej + log(uninfoCtr(i_index)./72);
end

percUninf = sum(uninfoCtr)/(numel(blocks).*n_trials).*100;
disp(['For this subject, ' num2str(sum(uninfoCtr)) ' uninformative trials were found.']);
disp([num2str(percUninf) ' percent of choice across all blocks uninformative.']);
disp(['All blocks (' num2str(numel(blocks)) ') simulated.']);

% Files to export to base workspace when called from another function
assignin('base', 'L_acc', L_acc);
assignin('base', 'L_rej', L_rej);
assignin('base','percUninf',percUninf);

objFn = ll_rej; % objective function is to minimize rejections and uninformative choice
end