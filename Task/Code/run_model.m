function [pcor,breg,lambda] = run_model(blck,sd_inf,nsim)
%  RUN_MODEL  Run model on ACTOBS experiment block
%
%  Usage: [pcor,breg,lambda] = RUN_MODEL(blck,sd_inf,nsim)
%
%  where blck   - block substructure
%        sd_inf - inference noise s.d.
%        nsim   - number of simulations (default:1e3)
%
%  The function returns the predicted proportion of correct responses pcor,
%  the maximum-likelihood parameter estimates breg of a logistic regression
%  with the following 3 regressors: 1/ the sequence evidence (in logLR units),
%  2/ the previous response (scaled as +1:left or -1:right) and 3/ a constant
%  term capturing response bias, and the best lambda parameter w.r.t. block-
%  wise accuracy.
%
%  Valentin Wyart <valentin.wyart@ens.fr> - 09/2015

% check input arguments
if nargin < 3
    nsim = 1e3;
end
if nargin < 2
    error('Missing input arguments!');
end

% get block-wise information
seqdir = blck.seqdir; % generative sequence laterality (1:left or 2:right)
seqllr = blck.seqllr; % sequence evidence (favoring left if >0 or right if <0)
seqlen = blck.seqlen; % sequence length
nseq = length(seqdir); % number of sequences

% simulate noisy inference
xs = seqllr; % signal
xn = sd_inf*sqrt(seqlen-1); % noise s.d.
x  = bsxfun(@plus,xs,bsxfun(@times,xn,randn(nsim,nseq))); % noisy signal

% find best lambda w.r.t. block-wise accuracy
lambda_vec = 0:0.01:1;
pcor = nan(size(lambda_vec));
for i = 1:length(lambda_vec)
    lambda = lambda_vec(i);
    % simulate model
    y = zeros(nsim,1); % accumulated evidence
    r = zeros(nsim,nseq); % responses (1:left or 2:right)
    for iseq = 1:nseq
        y = y*lambda+x(:,iseq);
        r(:,iseq) = 1+(y < 0);
    end
    % get accuracy
    pcor(i) = mean(mean(bsxfun(@eq,r,seqdir),2),1);
end
lambda = fminbnd(@(l)-interp1(lambda_vec,pcor,l,'spline'),0,1);

% re-simulate model with best lambda
y = zeros(nsim,1);
r = zeros(nsim,nseq+1);
r(:,1) = rand(nsim,1) < 0.5; % random 1st response
for iseq = 1:nseq
    y = y*lambda+x(:,iseq);
    r(:,iseq+1) = 1+(y < 0);
end
pcor = mean(mean(bsxfun(@eq,r(:,2:end),seqdir),2),1);

% perform logistic regression
xreg_all = [];
yreg_all = [];
for isim = 1:nsim
    xreg_sim = cat(2,seqllr',3-2*r(isim,1:nseq)');
    yreg_sim = r(isim,2:end)' == 1;
    % append current simulation
    xreg_all = cat(1,xreg_all,xreg_sim);
    yreg_all = cat(1,yreg_all,yreg_sim);
end
% add intercept term
xreg_all = cat(2,xreg_all,ones(size(xreg_all,1),1));
breg = logreg(xreg_all,yreg_all,'probit');
    
end