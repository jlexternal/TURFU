function [blck] = gen_blck(cfg)
%  GEN_BLCK  Generate block of ACTOBS experiment
%
%  Usage: [blck] = GEN_BLCK(cfg)
%
%  where the configuration structure cfg can contain the following fields:
%    * condtn - task condition => 1:stable or 2:volatile or 3:practice
%    * epimap - episode mapping => 1:orange=good|left or 2:blue=good|left
%    * epidir - initial episode direction => 1:left or 2:right
%
%  Valentin Wyart <valentin.wyart@ens.fr> - 09/2015

if nargin < 1
    error('Missing configuration structure!');
end

% get configuration parameters
condtn = cfg.condtn; % task condition
epimap = cfg.epimap; % episode mapping
epidir = cfg.epidir; % initial episode direction

% set block parameters
if condtn == 1
    % stable contingencies
    davg  = 12;
    nepi  = 6;
    kappa = 0.5;
elseif condtn == 2
    % volatile contingencies
    davg  = 6;
    nepi  = 12;
    kappa = 0.5;
elseif condtn == 3
    % practice
    davg  = 9;
    nepi  = 6;
    kappa = 1.0;
end
dlim = [4,24]; % minimum/maximum episode length

% number of sequences to be generated
nseq = nepi*davg;

% create block structure
blck           = [];
blck.condtn    = condtn;
blck.epimap    = epimap;
blck.cfg.davg  = davg;
blck.cfg.dlim  = dlim;
blck.cfg.nepi  = nepi;
blck.cfg.nseq  = nseq;
blck.cfg.kappa = kappa;

% set function handles
shuffle = @(x)x(randperm(numel(x)));
rndvm   = @(n,t,k)mod(randraw('vonmises',[t*2,k],n)/2,pi);
pdfvm   = @(x,t,k)exp(k.*cos(2*(x-t)))./(pi*besseli(0,k));
logit   = @(p)log(p./(1-p));
sigmd   = @(x)1./(1+exp(-x));

% set internal parameters
tprec = 0.1; % precision on theta (%)
kprec = 0.1; % precision on kappa (%)
nshuf = 1e4; % number of shuffles

scondtn = {'stable','volatile','practice'};
fprintf('Generating %s block:\n',scondtn{condtn});

% build episode length
fprintf('  * building list of episode lengths...\n');
out = gen_draws(davg,dlim,nepi);
epilen = out.xs;

% build episode direction
blck.epinum = []; % episode number
blck.epipos = []; % sequence position wrt episode onset
blck.epineg = []; % sequence position wrt episode offset
blck.epidir = []; % episode direction (1:good=left 2:good=right)
for i = 1:length(epilen)
    blck.epinum = cat(2,blck.epinum,i*ones(1,epilen(i)));
    blck.epipos = cat(2,blck.epipos,1:+1:epilen(i));
    blck.epineg = cat(2,blck.epineg,epilen(i):-1:1);
    blck.epidir = cat(2,blck.epidir,(mod(i+epidir,2)+1)*ones(1,epilen(i)));
end

nseq = sum(epilen); % number of sequences
nrev = nepi-1; % number of reversals

% build sequence direction
blck.seqdir = blck.epidir;
% remove episode direction = sequence direction
blck = rmfield(blck,'epidir');

% build sequence length
fprintf('  * building list of sequence lengths...\n');
seqlen0 = [2,4,6,8];
nlen = length(seqlen0);
nrep = floor((nepi-1)/nlen);
jepi = nepi-nrep*nlen+1:nepi;
while true
    retry = false;
    seqlen = zeros(size(blck.seqdir));
    % fill around reversals
    irev = [];
    for iepi = jepi
        irev = cat(1,irev,find(blck.epinum == iepi,1)+[0:3]);
    end
    nrev = size(irev,1);
    lrev = [];
    for i = 1:4
        while true
            lrev_tmp = [lrev,shuffle(kron(seqlen0(:),ones(nrep,1)))];
            isok = true;
            for j = 1:nrev
                isok = isok && ~HasConsecutiveValues(lrev_tmp(j,:),3);
            end
            if isok
                lrev = lrev_tmp;
                break
            end
        end
    end
    seqlen(irev(:)) = lrev(:);
    % fill gaps
    igap = find(seqlen == 0);
    ngap = length(igap);
    ndup = ceil(ngap/nlen);
    ntry = 0;
    while true
        ntry = ntry+1;
        if ntry > 1000
            retry = true;
            break
        end
        seqlen_tmp = seqlen;
        seqlen_tmp(igap) = randsample(kron(seqlen0,ones(1,ndup)),ngap,false);
        if ~HasConsecutiveValues(seqlen_tmp,3,[0])
            seqlen = seqlen_tmp;
            break
        end
    end
    if retry
        continue
    else
        break
    end
end
blck.seqlen = seqlen;

% build sequences
xt_all = cell(nlen,1);
llr_all = cell(nlen,1);
iseq_all = cell(nlen,1);
for ilen = 1:nlen
    fprintf('  * building sequences with %d samples...\n',seqlen0(ilen));
    nsamp = seqlen0(ilen);
    iseq = find(blck.seqlen == seqlen0(ilen));
    nseq = length(iseq);
    % generate samples
    while true
        xt = rndvm([nseq,nsamp],0,kappa);
        phat = fminsearch(@(p)-sum(log(pdfvm(xt(:),p(1),p(2)))),[0,kappa], ...
            optimset('Display','off','MaxFunEvals',1e4,'TolX',1e-6));
        that = phat(1);
        khat = phat(2);
        if abs(that)/pi < tprec/2 && abs(khat-kappa)/kappa < kprec/2
            break
        end
    end
    % group samples into sequences
    xt_perm = nan(nseq,nsamp,nshuf);
    llr_perm = nan(nseq,nshuf);
    for ishuf = 1:nshuf
        xt = reshape(xt(randperm(nseq*nsamp)),[nseq,nsamp]);
        llr = sum(2*kappa*cos(2*xt),2);
        [~,isort] = sort(llr);
        xt_perm(:,:,ishuf) = xt(isort,:);
        llr_perm(:,ishuf) = llr(isort);
    end
    llr_diff = bsxfun(@minus,llr_perm,mean(llr_perm,2));
    [~,imin] = min(sum(abs(llr_diff),1));
    xt = xt_perm(:,:,imin);
    llr = llr_perm(:,imin);
    clear xt_perm llr_perm llr_diff
    xt_all{ilen} = xt;
    llr_all{ilen} = llr;
    iseq_all{ilen} = iseq;
end

% shuffle sequences
fprintf('  * building list of sequences...\n');
ishuf_perm = cell(nlen,nshuf);
nseq = blck.cfg.nseq;
ps_perm = nan(nshuf,nseq);
lls = 2*(blck.seqdir == 1)-1;
for ishuf = 1:nshuf
    % shuffle sequences
    llr = zeros(1,nseq);
    for ilen = 1:nlen
        ishuf_perm{ilen,ishuf} = randperm(length(iseq_all{ilen}));
        llr(iseq_all{ilen}) = llr_all{ilen}(ishuf_perm{ilen,ishuf});
    end
    llr = llr.*lls;
    % simulate model 
    pr = out.pb; % p(reversal)
    ps = zeros(1,nseq);
    ps(1) = 0.5;
    for i = 1:nseq
        ps(i) = (1-pr)*ps(i)+pr*(1-ps(i));
        ps(i) = sigmd(logit(ps(i))+llr(i));
        if i < nseq
            ps(i+1) = ps(i);
        end
    end
    ps_perm(ishuf,:) = ps;
end
qs_perm = logit(ps_perm);
zs_perm = zscore(qs_perm,[],1);
[~,imin] = min(sum(zs_perm.^2,2));

% add sequences to experiment
seqtlt = cell(size(blck.seqdir)); % list of sample tilts w.r.t. category center
seqllr = nan(size(blck.seqdir)); % evidence in favor of left response
ishuf = ishuf_perm(:,imin);
for ilen = 1:nlen
    iseq = iseq_all{ilen};
    xt = xt_all{ilen}(ishuf{ilen},:);
    llr = llr_all{ilen}(ishuf{ilen});
    for i = 1:length(iseq)
        seqtlt{iseq(i)} = xt(i,:);
        if blck.seqdir(iseq(i)) == 1
            seqllr(iseq(i)) = +llr(i);
        else
            seqllr(iseq(i)) = -llr(i);
        end
    end
end
blck.seqtlt = seqtlt;
blck.seqllr = seqllr;

% add fields that need to be filled online
blck.seqcat = zeros(size(blck.seqdir)); % sequence category
blck.seqang = cell(size(blck.seqtlt)); % list of sample angles

% shuffle category angles
fprintf('  * building list of category mappings...\n');
catmap = [];
for i = 1:ceil(nseq/8)
    while true
        catmap_tmp = cat(2,catmap,randperm(8));
        if ...
                ~HasConsecutiveValues(catmap_tmp,3) && ...
                ~HasConsecutiveValues(ceil(catmap_tmp/4),4)
            break
        end
    end
    catmap = catmap_tmp;
end
catmap = catmap(1:nseq);
blck.catmap = catmap;

fprintf('Done generating block.\n\n');

end