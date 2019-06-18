function [MU,A,b,V,n] = updateclusters(SS0,SS1,SS2,pr)
%__________________________________________________________________________
% anatomy.segment.gmm.lib.updateclusters
%--------------------------------------------------------------------------
% FORMAT [MU,A,b,V,n] = gmm.lib.updateclusters(SS0,SS1,SS2,{MU0,b0,V0,n0})
% SS0 - 0th order sufficient statistics (sum Z_i)
% SS1 - 1st order sufficient statistics (sum Z_i * X_i)
% SS2 - 2nd order sufficient statistics (sum Z_i * (X_i * X_i'))
% pr  - List of prior Gauss-Wishart parameters.
%
% Compute posterior GMM parameters from suff stats.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

import anatomy.math.matrix.*       % posdef.inv

if nargin<4, pr =[]; end

K  = numel(SS0);
MU0 = [];
b0  = [];
V0  = [];
n0  = [];
if numel(pr) >= 1
    MU0 = pr{1};
    if numel(pr) >= 2
        b0 = pr{2};
        if numel(pr) >= 3
            V0 = pr{3};
            if numel(pr) >= 4
                n0 = pr{4};
            end
        end
    end
end

%--------------------------------------------------------------------------
% Use double precision
MU0 = double(MU0);
b0  = double(b0);
V0  = double(V0);
n0  = double(n0);

% -------------------------------------------------------------------------
% Mean
if sum(b0) == 0
    % ---------------------------------------------------------------------
    % Without prior
    b  = [];
    MU = bsxfun(@rdivide, SS1, SS0);
else
    % ---------------------------------------------------------------------
    % With prior
    b  = b0 + SS0;
    MU = bsxfun(@rdivide, SS1 + bsxfun(@times,b0,MU0), b);
end

% -------------------------------------------------------------------------
% Scale/Precision
if sum(n0) == 0
    % ---------------------------------------------------------------------
    % Without prior
    n   = [];
    for k=1:K
        SS2(:,:,k) = SS2(:,:,k) / SS0(k) - (MU(:,k) * MU(:,k).');
    end
else
    % ---------------------------------------------------------------------
    % With prior
    n = n0 + SS0;
    for k=1:K
        SS2(:,:,k) = SS2(:,:,k) +   b0(k) * MU0(:,k) * MU0(:,k).' ...
                                -    b(k) * MU(:,k)  * MU(:,k).' ...
                                + posdef.inv(V0(:,:,k));
    end
end
V = SS2;
for k=1:K
    V(:,:,k) = posdef.inv(V(:,:,k));
end
if sum(n) > 0
    A = bsxfun(@times, V, reshape(n, [1 1 K]));
else
    A = V;
    V = [];
end