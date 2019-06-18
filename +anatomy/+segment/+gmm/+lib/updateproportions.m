function [PI,logPI,a] = updateproportions(SS0, a0)
%__________________________________________________________________________
% anatomy.segment.gmm.lib.updateproportions
%--------------------------------------------------------------------------
% FORMAT [PI,logPI,a,ll] = gmm.lib.updateproportions(SS0, a0)
%
% SS0 - 1xK 0th order sufficient statistics (sum of responsibilities)
% a0  - 1xK Dirichlet prior (can be 0)
%
% PI    - 1xK Cluster proportion posterior expected value
% logPI - 1xK ln(PI) or E[ln(PI)] (if Bayesian)
% a     - Dirichlet posterior (if Bayesian)
% ll    - 1x1 Lower bound: E[ln p(PI|a)] - E[ln q(PI)]
%
% Bayesian or ML update of cluster proportions.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

a = a0 + SS0;
if sum(a0(:))
% Bayesian
    % expected values
    logPI = psi(a) - psi(sum(a));
    PI    = a ./ sum(a(:));
else
% Maximum Likelihood
    a     = max(a, eps);
    PI    = a ./ sum(a(:));
    logPI = log(PI);
end