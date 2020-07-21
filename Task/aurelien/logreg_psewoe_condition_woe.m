function [pparam,lps,llh,aic,bic] = logreg_psewoe_condition_woe(x,y,flink,ifilt)

if nargin < 4
    flink = 'logit';
end

if nargin < 3
    error('Missing input arguments!');
end

% ensure binary output array
y = logical(y);

% setup gradient-descent function
options = optimset('fmincon');
options = optimset(options,'Algorithm','interior-point','Display','notify', ...
    'FunValCheck','on','TolX',1e-20,'MaxFunEvals',1e6);

% run maximum-likelihood estimation

    [pval,fval] = fmincon(@(p)-get_llh(p), ...
        [1,1,1,0,0],[],[],[],[],[0,0,0,0,0],[6,6,6,0.5,0.5],[],options);
pparam = [pval(1),pval(3);pval(2),pval(2)];
lps    = [pval(4),pval(5)];

nobs = size(x,1);
nfit = 2;
llh = get_llh(pval);
aic = -2*llh+2*nfit+2*nfit*(nfit+1)/(nobs-nfit+1);
bic = -2*llh+nfit*log(nobs);

    function [llh] = get_llh(p)
        dv(ifilt == 1) = x(ifilt == 1).*p(2) + p(1);
        dv(ifilt == 2) = x(ifilt == 2).*p(2) + p(3);
        
        switch flink
            case 'logit'
                plh = 1./(1+exp(-dv));
            case 'probit'
                plh = normcdf(dv);
        end
        plh(ifilt == 1) = p(4) + (1-p(4))*plh(ifilt == 1);
        plh(ifilt == 2) = p(5) + (1-p(5))*plh(ifilt == 2);
        plh(y==0) = 1-plh(y==0);
        llh = sum(log(max(plh,realmin)));
    end

end