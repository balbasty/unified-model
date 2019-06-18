function X = missing(X, Z, cluster, codes)
%__________________________________________________________________________
% anatomy.segment.gmm.lib.missing
%--------------------------------------------------------------------------
% FORMAT X = gmm.lib.missing(X, Z, {MU,A}, {C,L})
% X  - NxP   observations
% Z  - NxK   responsibilities
% MU - PxK   (expected) means
% A  - PxPxK (expected) precision matrices
% C  - Nx1   "missing value" code image
% L  -       list of existing codes
%
% X - NxP    observations with inferred values
%
% Compute the mean expected value of missing voxels.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

import anatomy.segment.gmm.lib.code2bin

MU = [];
A  = [];
C  = [];
L  = [];

%--------------------------------------------------------------------------
% Read input arguments
if ~iscell(cluster)
    MU = cluster;
else
    if numel(cluster) >= 1
        MU = cluster{1};
        if numel(cluster) >= 2
            A = cluster{2};
        end
    end
end
if nargin >= 4
    if ~iscell(codes)
        C = codes;
    else
        if numel(codes) >= 1
            C = codes{1};
            if numel(codes) >= 2
                L = codes{2};
            end
        end
    end
    if isempty(L)
        L = unique(C);
    end
end

%--------------------------------------------------------------------------
% Dimensions
P = size(X, 2);
K = size(Z, 2);
if isempty(L)
    L = 2^P - 1; % None missing
end

% -------------------------------------------------------------------------
% For each missing combination
for i=1:numel(L)
    
    % ---------------------------------------------------------------------
    % Get code / mask / missing modalities
    c        = L(i);
    observed = code2bin(c, P);
    missing  = ~observed;
    Pm       = sum(missing);
    if Pm == 0, continue; end
    msk      = (C == c);
    Nm       = sum(msk);
    if Nm == 0, continue; end
    
    % ---------------------------------------------------------------------
    % Initialise
    X(msk,missing) = 0;

    % ---------------------------------------------------------------------
    % Compute posterior mean (expected value)
    % 1) t = sum_k {z * ( mu[m] + A[m]/A[m,o]*(mu[o]-g) ) }
    for k=1:K
        X1k = zeros(1, 'like', X);
        X1k = bsxfun(@plus,X1k,MU(missing,k).');
        X1k = bsxfun(@plus,X1k,bsxfun(@minus, MU(observed,k).', X(msk,observed)) * (A(observed,missing,k) / A(missing,missing,k)));
        X(msk,missing) = X(msk,missing) + bsxfun(@times, X1k, Z(msk,k));
    end
                    
end