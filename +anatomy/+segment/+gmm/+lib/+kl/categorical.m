function klZ = categorical(Z, W, logPI)
%__________________________________________________________________________
% anatomy.segment.gmm.lib.kl.categorical
%--------------------------------------------------------------------------
% klZ = math.gmm.lib.kl.categorical(Z, W, logPI)
% 
% KL divergence between two Categorical distributions
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging


% Initialise
klZ = zeros(1, 'like', Z);

% E[ln p(Z|PI)] (prior ~ responsibilities)
klZ = klZ + sum(sum(bsxfun(@times,Z,logPI), 2) .* W, 1);

% -E[ln q(Z)] (posterior ~ responsibilities))
klZ = klZ - sum(sum(Z .* log(max(Z,eps)), 2) .* W, 1);
