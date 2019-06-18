function lg = loggamma(a, p)
%__________________________________________________________________________
% anatomy.math.probability.loggamma
%--------------------------------------------------------------------------
% lg = loggamma(a, p)
%  
% Log of multivariate gamma function of order p
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging
    if nargin < 2
        p = 1;
    end
    lg = (p*(p-1)/4)*log(pi);
    for i=1:p
        lg = lg + gammaln(a + (1-p)/2);
    end
end