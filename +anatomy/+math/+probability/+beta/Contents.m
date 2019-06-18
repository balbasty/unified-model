%__________________________________________________________________________
% anatomy.math.probability.beta
%--------------------------------------------------------------------------
% Characteristic functions of the Beta distribution:
%
%   [pdf]       Probability density function
%   [ll/logpdf] Log-probability density function
%   [kl]        Kullback-Leibler divergence
%   [up/update] Conjugate (or Bayesian) update
%   [E]         Expected value (E[x])
%   [Elog]      Expected log value (E[ln x])
%   [V]         Variance (V[x])
%   [Vlog]      Variance of the log value (V[ln x])
%
%--------------------------------------------------------------------------
% Relationship to other distributions
% -----------------------------------
%
% The Beta distribution is a conjugate prior for the proportion parameter 
% of the Bernoulli, Binomial, negative Binomial and Geometric 
% distributions. 
%
% The beta prameter is a special case of the Dirichlet distribution when
% the number of classes is K = 2, i.e., the underlying trial is binary.
%
%--------------------------------------------------------------------------
% General distribution
% --------------------
%
% The Beta distribution is parameterised by two strictly positive shape 
% parameters, a and b, and is defined over [0,1]. Parameters a and b can be
% seen as concentration parameters, i.e., they reflect both class
% proportions and their variance.
% The expected value of the distribution is E[x] = a/(a+b)
% The variance of the distribution is       V[x] = ab/[(a+b)^2*(a+b+1)]
%
% FORMAT pdf = beta.pdf(x, a, b)
% FORMAT ll  = beta.logpdf(x, a, b)
%   >> (Log) Probability density function.
%
% FORMAT e  = beta.E(a, b)
% FORMAT el = beta.Elog(a, b)
% FORMAT v  = beta..V(a, b)
% FORMAT vl = beta.Vlog(a, b)
%   >> Mean and variance
%
% FORMAT kl  = beta.kl(a1, b1, a0, b0)
%   >> Kullback-Leibler divergence from B0 to B1 = KL(B1||B0).
%
%--------------------------------------------------------------------------
% Bernoulli probability conjugate
% -------------------------------
%
% To make the distribution easier to manipulate when used as a conjugate
% prior for a Bernoulli distribution, we reparameterise it by setting:
% - a = np
% - b = n(1-p)
% The expected value of the distribution becomes E[x] = p
% The variance of the distribution becomes       V[x] = p(1-p)/(n+1)
% Then, the Beta distribution is parameterised by a mean probability  
% parameter (p) and a degrees of freedom (n).
%
% FORMAT pdf = beta.pdf(x, p, n, 'ber')
% FORMAT ll  = beta.logpdf(x, p, n, 'ber')
%   >> (Log) Probability density function.
%
% FORMAT e  = beta.E(p, n, 'ber')
% FORMAT el = beta.Elog(p, n, 'ber')
% FORMAT v  = beta.V(p, n, 'ber')
% FORMAT vl = beta.Vlog(p, n, 'ber')
%   >> Mean and variance.
%
% FORMAT kl  = beta.kl(p1, n1, p0, n0, 'ber')
%   >> Kullback-Leibler divergence from B0 to B1 = KL(B1||B0).
%
% FORMAT [p1, n1] = beta.update(p,  n,  p0, n0)
% FORMAT [p1, n1] = beta.update(s0, s1, p0, n0, 'ber')
%   >> Posterior parameters of the Beta distribution.
%
%--------------------------------------------------------------------------
% Binomial probability conjugate
% ------------------------------
%
% To make the distribution easier to manipulate when used as a conjugate
% prior for a Binomial distribution, we reparameterise it by setting:
% - a = knp
% - b = kn(1-p)
% The expected value of the distribution becomes E[x] = p
% The variance of the distribution becomes       V[x] = p(1-p)/(kn+1)
% Then, the Beta distribution is parameterised by a mean probability  
% parameter (p), a degrees of freedom (n) and a number of trials (k).
%
% FORMAT pdf = beta.pdf(x, p, n, k, 'bin')
% FORMAT ll  = beta.logpdf(x, p, n, k, 'bin')
%   >> (Log) Probability density function.
%
% FORMAT e  = beta.E(p, n, k, 'bin')
% FORMAT el = beta.Elog(p, n, k, 'bin')
% FORMAT v  = beta.V(p, n, k, 'bin')
% FORMAT vl = beta.Vlog(p, n, k, 'bin')
%   >> Mean and variance.
%
% FORMAT kl  = beta.kl(p1, n1, p0, n0, k, 'bin')
%   >> Kullback-Leibler divergence from B0 to B1 = KL(B1||B0).
%
% FORMAT [p1, n1] = beta.update(p,  n,  p0, n0)
% FORMAT [p1, n1] = beta.update(s0, s1, p0, n0, k, 'bin')
%   >> Posterior parameters of the Beta distribution.
%
%--------------------------------------------------------------------------
% Negative-Binomial probability conjugate
% ---------------------------------------
%
% To make the distribution easier to manipulate when used as a conjugate
% prior for a negative Binomial distribution, we reparameterise it by 
% setting:
% - a = rnp/(1-p)
% - b = rn
% The expected value of the distribution becomes E[x] = p
% The variance of the distribution becomes       V[x] = p(1-p)^2/(rn+1-p)
% Then, the Beta distribution is parameterised by a mean probability  
% parameter (p), a degrees of freedom (n) and a number of failures (r).
%
% FORMAT pdf = beta.pdf(x, p, n, r, 'nbin')
% FORMAT ll  = beta.logpdf(x, p, n, r, 'nbin')
%   >> (Log) Probability density function.
%
% FORMAT e  = beta.E(p, n, r, 'nbin')
% FORMAT el = beta.Elog(p, n, r, 'nbin')
% FORMAT v  = beta.V(p, n, r, 'nbin')
% FORMAT vl = beta.Vlog(p, n, r, 'nbin')
%   >> Mean and variance.
%
% FORMAT kl  = beta.kl(p1, n1, p0, n0, r, 'nbin')
%   >> Kullback-Leibler divergence from B0 to B1 = KL(B1||B0).
%
% FORMAT [p1, n1] = beta.update(p,  n,  p0, n0)
% FORMAT [p1, n1] = beta.update(s0, s1, p0, n0, r, 'nbin')
%   >> Posterior parameters of the Beta distribution.
%
%--------------------------------------------------------------------------
% Geometric probability conjugate
% -------------------------------
%
% To make the distribution easier to manipulate when used as a conjugate
% prior for a Geometric distribution, we reparameterise it by 
% setting:
% - a = n
% - b = n(1-p)/p
% The expected value of the distribution becomes E[x] = p
% The variance of the distribution becomes       V[x] = (1-p)/(n+p)
% Then, the Beta distribution is parameterised by a mean probability  
% parameter (p) and a degrees of freedom (n).
%
% FORMAT pdf = beta.pdf(x, p, n, 'geom')
% FORMAT ll  = beta.logpdf(x, p, n, 'geom')
%   >> (Log) Probability density function.
%
% FORMAT e  = beta.E(p, n, 'geom')
% FORMAT el = beta.Elog(p, n, 'geom')
% FORMAT v  = beta.V(p, n, 'geom')
% FORMAT vl = beta.Vlog(p, n, 'geom')
%   >> Mean and variance.
%
% FORMAT kl  = beta.kl(p1, n1, p0, n0, 'geom')
%   >> Kullback-Leibler divergence from B0 to B1 = KL(B1||B0).
%
% FORMAT [p1, n1] = beta.update(p,  n,  p0, n0)
% FORMAT [p1, n1] = beta.update(s0, s1, p0, n0, 'geom')
% FORMAT [p1, n1] = beta.update(s0, s1, p0, n0, 'fgeom')
%   >> Posterior parameters of the Beta distribution.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging