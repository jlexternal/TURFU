function [pres] = gen_pres(seqcat)

% check input arguments
if nargin < 1
    error('Missing sequence category!');
end

addpath('./Toolboxes/Rand');

% setup random number generator
SetupRand;

% set parameters
kappa = 0.5; % kappa
tprec = 0.1; % precision on theta (%)
kprec = 0.1; % precision on kappa (%)
nshuf = 1e4; % number of shuffles

% set function handles
shuffle = @(x)x(randperm(numel(x)));
rndvm   = @(n,t,k)mod(randraw('vonmises',[t*2,k],n)/2,pi);
pdfvm   = @(x,t,k)exp(k.*cos(2*(x-t)))./(pi*besseli(0,k));
logit   = @(p)log(p./(1-p));
sigmd   = @(x)1./(1+exp(-x));

seqlen0 = [2,4,6,8]; % list of sequence lengths
nseq = 4; % number of sequences per length

nlen = length(seqlen0);
ntot = nseq*nlen;
if mod(ntot,8) > 0
    error('Total number of sequences must be divisible by 8!');
end

% generate sequences
pres        = [];
pres.seqlen = kron(seqlen0,ones(1,nseq));
pres.seqtlt = {};
pres.seqllr = [];
for ilen = 1:nlen
    nsamp = seqlen0(ilen);
    % generate samples
    xt = rndvm([nseq*nsamp,nshuf],0,kappa);
    xt = reshape(xt,[nseq,nsamp,nshuf]);
    % pick most average draw w.r.t. sorted sequence-wise evidence
    llr = squeeze(sum(2*kappa*cos(2*xt),2));
    llr = sort(llr,1);
    [~,ishuf] = min(sum(bsxfun(@minus,llr,mean(llr,2)).^2,1));
    xt = xt(:,:,ishuf);
    llr = llr(:,ishuf);
    for iseq = 1:nseq
        pres.seqtlt{1,end+1} = xt(iseq,:);
        pres.seqllr(1,end+1) = llr(iseq);
    end
end

% shuffle sequences
while true
    iperm = randperm(ntot);
    if ...
            ~HasConsecutiveValues(pres.seqlen(iperm),2) && ...
            ~HasConsecutiveValues(mod(iperm,nseq),2)
        break
    end
end
pres.seqlen = pres.seqlen(iperm);
pres.seqtlt = pres.seqtlt(iperm);
pres.seqllr = pres.seqllr(iperm);
pres.seqcat = seqcat(ones(size(pres.seqlen)));

% shuffle category mappings
pres.catmap = kron(1:8,ones(1,ntot/8));
while true
    pres.catmap = shuffle(pres.catmap);
    if ...
            ~HasConsecutiveValues(pres.catmap,2) && ...
            ~HasConsecutiveValues(ceil(pres.catmap/4),4)
        break
    end
end
    
end