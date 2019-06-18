%__________________________________________________________________________
% anatomy.math.probability.normal
%--------------------------------------------------------------------------
% Characteristic functions of the (uni/multivariate) Normal distribution:
%   [pdf]       Probability density function
%   [ll/logpdf] Log-probability density function
%   [kl]        Kullback-Leibler divergence
%   [up/update] Conjugate (or Bayesian) update
%
%--------------------------------------------------------------------------
% General distribution
% --------------------
%
% The Normal distribution is parameterised by a mean parameter (mu) and
% a (co)variance (sigma) or precision (lambda) parameter. 
% For multivariate distributions, sigma is a KxK covariance matrix. 
% In the univariate case, it reduces to a scalar value, the variance.
%
% FORMAT pdf = normal.pdf(x, mu, sigma)
% FORMAT pdf = normal.pdf(x, mu, lambda, 'precision')
% FORMAT ll  = normal.logpdf(x, mu, sigma)
% FORMAT ll  = normal.logpdf(x, mu, lambda, 'precision')
%   >> (Log) Probability density function.
%
% FORMAT kl  = normal.kl(mu1, sigma1,  mu0, sigma0)
% FORMAT kl  = normal.kl(mu1, lambda1, mu0, lambda0, 'precision')
%   >> Kullback-Leibler divergence from N0 to N1 = KL(N1||N0)
%
%--------------------------------------------------------------------------
% Normal mean conjugate
% ---------------------
%
% The Normal distribution can be used as a conjugate prior for the mean
% parameter of another Normal distribution with known covariance.
% It is then parameterised by an expected mean (mu), a degrees of freedom 
% (n) and a known covariance (sigma) or precision (lambda). 
%
% FORMAT pdf = normal.pdf(x, mu, n, sigma)
% FORMAT pdf = normal.pdf(x, mu, n, lambda, 'precision')
% FORMAT ll  = normal.logpdf(x, mu, n, sigma)
% FORMAT ll  = normal.logpdf(x, mu, n, lambda, 'precision')
%   >> (Log) Probability density function.
%
% FORMAT kl  = normal.kl(mu1, n1, mu0, n0, sigma)
% FORMAT kl  = normal.kl(mu1, n1, mu0, n0, lambda, 'precision')
%   >> Kullback-Leibler divergence from N0 to N1 = KL(N1||N0)
%
% FORMAT [mu1, n1] = normal.update(mu, n, mu0, n0)
%   >> Posterior parameters of the Normal distribution.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging