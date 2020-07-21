function [b] = logreg(x,y,flink)
%  LOGREG  Perform binary logistic regression
%
%  Usage: [b] = LOGREG(x,y,flink)
%
%  where x     - design matrix (input)
%        y     - binary responses (output)
%        flink - link function ('logit' or 'probit')
%
%  The function returns maximum-likelihood parameter estimates b obtained via
%  gradient descent using fminsearch. This function serves as a minimalist
%  alternative to glmfit where the Statistics Toolbox is not available.
%
%  Valentin Wyart <valentin.wyart@ens.fr> - 09/2015

% check input arguments
if nargin < 3
    error('Missing input arguments!');
end
if ~ismember(flink,{'logit','probit'})
    error('Invalid link function %s!',flink);
end
if size(x,1) ~= size(y,1)
    error('Mismatching input/output arrays!');
end

% ensure binary output array
y = logical(y);

% get maximum-likelihood parameter estimates
n = size(x,2);
b = fminsearch(@(pp)-get_llh(pp),zeros(n,1));

    function [llh] = get_llh(p)
        % get logistic regression model log-likelihood
        dv = (x*p).*(y*2-1);
        switch flink
            case 'logit'
                plh = 1./(1+exp(-dv));
            case 'probit'
                plh = normcdf(dv);
        end
        llh = sum(log(max(plh,eps)));
    end

end