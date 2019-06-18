function [mu1, n1] = update(mu, n, mu0, n0)
%__________________________________________________________________________
% anatomy.math.probability.normal.update
%--------------------------------------------------------------------------
% FORMAT [mu1, n1] = update(mu, n, mu0, n0)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    n1 = n + n0;
    mu1 = (n0 * mu0 + n * mu)/n1;

end