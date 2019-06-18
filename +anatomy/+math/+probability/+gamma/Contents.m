%__________________________________________________________________________
% anatomy.math.probability.gamma
%--------------------------------------------------------------------------
% Characteristic functions of the Gamma distribution:
%   [pdf]       Probability density function
%   [ll/logpdf] Log-probability density function
%   [kl]        Kullback-Leibler divergence
%   [up/update] Conjugate (or Bayesian) update
%   [E]         Expected value (E[x])
%   [Elog]      Expected log (E[ln x])
%   [V]         Variance (V[x])
%   [Vlog]      Variance of the log (V[ln x])
%
% The Gamma distribution is a conjugate prior for a Normal precision (or
% precision magnitude) with known mean, for a Gamma rate with known shape 
% or in general for any rate parameter of an Exponential family
% distribution.
%
%--------------------------------------------------------------------------
% General distribution
% --------------------
%
% The Gamma distribution is parameterised by a shape parameter (alpha) and
% a rate parameter (beta).
%
% FORMAT pdf = gamma.pdf(x, alpha, beta)
% FORMAT ll  = gamma.logpdf(x, alpha, alpha)
%   >> (Log) Probability density function.
%
% FORMAT e  = gamma.E(alpha, beta)
% FORMAT el = gamma.Elog(alpha, beta)
% FORMAT v  = gamma.V(alpha, beta)
% FORMAT vl = gamma.Vlog(alpha, beta)
%   >> Mean and variance
%
% FORMAT kl  = gamma.kl(alpha1, beta1, alpha0, beta0)
%   >> Kullback-Leibler divergence from G0 to G1 = KL(G1||G0).
%
%--------------------------------------------------------------------------
% Normal precision conjugate
% --------------------------
%
% The Gamma distribution is parameterised by a mean precision parameter 
% (lambda) and a degrees of freedom (n). It can be a precision *magnitude*
% parameter of K > 1.
%
% FORMAT pdf = gamma.pdf(x, lambda, n, K, ('normal'))
% FORMAT ll  = gamma.logpdf(x, lambda, n, K, ('normal'))
%   >> (Log) Probability density function.
%
% FORMAT e  = gamma.E(lambda, n, K, ('normal'))
% FORMAT el = gamma.Elog(lambda, n, K, ('normal'))
% FORMAT v  = gamma.V(lambda, n, K, ('normal'))
% FORMAT vl = gamma.Vlog(lambda, n, K, ('normal'))
%   >> Mean and variance.
%
% FORMAT kl  = gamma.kl(lam1, n1, lam0, n0, K, ('normal'))
%   >> Kullback-Leibler divergence from G0 to G1 = KL(G1||G0).
%
% FORMAT [lam1, n1] = gamma.update(lam, n,     lam0, n0, ('normal'))
% FORMAT [lam1, n1] = gamma.update(s0, s1, s2, lam0, n0, 
%                                  (mu=0), (Lam=eye), ('normal'))
%   >> Posterior parameters of the Gamma distribution.
%
%--------------------------------------------------------------------------
% Gamma rate conjugate
% --------------------
%
% The Gamma distribution is parameterised by a mean rate parameter 
% (beta), a degrees of freedom (n), and a known shape parameter (alpha).
%
% FORMAT pdf = gamma.pdf(x, beta, n, alpha, 'gamma')
% FORMAT ll  = gamma.logpdf(x, beta, n, alpha, 'gamma')
%   >> (Log) Probability density function.
%
% FORMAT e  = gamma.E(beta, n, alpha, 'gamma')
% FORMAT el = gamma.Elog(beta, n, alpha, 'gamma')
% FORMAT v  = gamma.V(beta, n, alpha, 'gamma')
% FORMAT vl = gamma.Vlog(beta, n, alpha, 'gamma')
%   >> Mean and variance.
%
% FORMAT kl  = gamma.kl(beta1, n1, beta0, n0, alpha, 'gamma')
%   >> Kullback-Leibler divergence from G0 to G1 = KL(G1||G0).
%
% FORMAT [beta1, n1] = gamma.update(beta, n, beta0, n0, 'gamma')
% FORMAT [beta1, n1] = gamma.update(s0, s1,  beta0, n0, alpha, 'gamma')
%   >> Posterior parameters of the Gamma distribution.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging