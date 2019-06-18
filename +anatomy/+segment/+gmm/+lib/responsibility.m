function Z = responsibility(logpX, logPI, varargin)
%__________________________________________________________________________
% anatomy.segment.gmm.lib.responsibility
%--------------------------------------------------------------------------
% FORMAT Z = gmm.lib.responsibility(logpX, logPI)
%
% Compute responsibilities.
% Responsibilities are the posterior expected value of class-indexing 
% vectors z_n.
% The posterior is computed as:
%   r_nk = exp(E[log Pi_k] + E[log p(x_n | Theta_k)]) / sum_k {r_nk}
%
% Extra terms can be added to the responsibilities prior to softmax by the
% varargin argument. These arguments need to be compatible (w.r.t. size) with 
% the following function call: bsxfun(@plus, Z, varargin{i}).
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

% Use double precision
logpX = double(logpX);
logPI = double(logPI);

% Omit NaN
logpX(isnan(logpX)) = 0; 

% Add terms: E[log Pi_k] + E[log p(X | Theta_k)]
Z = bsxfun(@plus, logpX, logPI);

for i=1:numel(varargin)
    Z = bsxfun(@plus, Z, double(varargin{i}));
end

% Exponentiate and normalise
Z = bsxfun(@minus, Z, max(Z, [], 2));
Z = exp(Z);
Z = bsxfun(@rdivide, Z, sum(Z, 2));