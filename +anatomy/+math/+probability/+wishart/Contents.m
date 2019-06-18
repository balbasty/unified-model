%__________________________________________________________________________
% anatomy.math.probability.wishart
%--------------------------------------------------------------------------
% Characteristic functions of the Wishart distribution:
%   [pdf]       Probability density function
%   [ll/logpdf] Log-probability density function
%   [kl]        Kullback-Leibler divergence
%   [up/update] Conjugate (or Bayesian) update
%   [E]         Expected value (E[X])
%   [Elogdet]   Expected log determinant (E[ln det X])
%   [V]         Variance (V[X])
%   [Vlogdet]   Variance of the log determinant (V[ln det X])
%
% The Wishart distribution is a conjugate prior for a multivariate Normal 
% precision matrix with known mean.
%
%--------------------------------------------------------------------------
% General distribution
% --------------------
%
% The Wishart distribution is parameterised by a scale matrix (V) and a
% degrees of freedom (n). It can be seen as the distribution of the sum of
% n independent centered multivariate Normal variables with precision 
% matrix V.
%
% FORMAT pdf  = wishart.pdf(X, V, n)
% FORMAT ll   = wishart.logpdf(X, V, n)
% FORMAT logZ = wishart.logconst(V, n)
%   >> (Log) Probability density function.
%
% FORMAT e    = wishart.E(V, n)
% FORMAT el   = wishart.Elogdet(V, n)
% FORMAT v    = wishart.V(V, n)
% FORMAT vl   = wishart.Vlogdet(V, n)
%   >> Mean and variance
%
% FORMAT kl  = wishart.kl(V1, n1, V0, n0)
%   >> Kullback-Leibler divergence from W0 to W1 = KL(W1||W0).
%
%--------------------------------------------------------------------------
% Normal precision matrix conjugate
% ---------------------------------
%
% The Wishart distribution is parameterised by a mean precision parameter 
% (Lambda) and a degrees of freedom (n).
%
% FORMAT pdf = wishart.pdf(X, Lambda, n, 'normal')
% FORMAT ll  = wishart.logpdf(X, Lambda, n, 'normal')
%   >> (Log) Probability density function.
%
% FORMAT e  = wishart.E(Lambda, n, 'normal')
% FORMAT el = wishart.Elogdet(Lambda, n, 'normal')
% FORMAT v  = wishart.V(Lambda, n, 'normal')
% FORMAT vl = wishart.Vlogdet(Lambda, n, 'normal')
%   >> Mean and variance.
%
% FORMAT kl  = wishart.kl(Lam1, n1, Lam0, n0, 'normal')
%   >> Kullback-Leibler divergence from W0 to W1 = KL(W1||W0).
%
% FORMAT [lam1, n1] = wishart.update(Lam, n,     Lam0, n0)
% FORMAT [lam1, n1] = wishart.update(s0, s1, s2, Lam0, n0, (mu=0))
%   >> Posterior parameters of the Wishart distribution.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging