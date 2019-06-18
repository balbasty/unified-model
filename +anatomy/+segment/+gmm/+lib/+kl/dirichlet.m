function klP = dirichlet(a, a0)
%__________________________________________________________________________
% anatomy.segment.gmm.lib.kl.dirichlet
%--------------------------------------------------------------------------
% klP = gmm.lib.kl.Dirichlet(a, a0)
% 
% KL divergence between two Dirichlet distributions
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

klP = zeros(1, 'like', a);
K   = numel(a);
if sum(a0) > 0
    % prior
    klP = gammaln(sum(a0)) - sum(gammaln(a0));
    klP = klP + sum((a0-1) .* (psi(a) - K*psi(sum(a))));
    % posterior
    klP = klP - gammaln(sum(a)) - sum(gammaln(a));
    klP = klP - sum((a-1) .* (psi(a) - K*psi(sum(a))));
end
