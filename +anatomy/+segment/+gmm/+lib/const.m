function c = const(mean,prec,L)
%__________________________________________________________________________
% anatomy.segment.gmm.lib.const
%--------------------------------------------------------------------------
% FORMAT c = gmm.lib.const({MU,b}, {V,n})
% FORMAT c = gmm.lib.const({MU,b}, {A})
% FORMAT c = gmm.lib.const({MU},   {A})
% FORMAT c = gmm.lib.const(..., L)
% MU - (Expected) mean
% b  - Mean df (if isempty or 0 -> no Bayesian prior)
% V  - Scale matrix     (if not n isempty or 0) 
% A  - Precision matrix (if n isempty or 0)
% n  - Precision df (if isempty or 0 -> no Bayesian prior)
%
% L - If provided, compute one term for each combination of 
%             missing data
%
% Compute the constant term (w.r.t. voxels) of each Gaussian 
% (expected) log-likelihood.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

import anatomy.segment.gmm.lib.code2bin
import anatomy.math.matrix.*           % posdef.logdet
import anatomy.math.probability.*      % wishart.Elogdet

MU = [];
b  = [];
V  = []; % It can actually be A (when n == 0)
n  = [];

%--------------------------------------------------------------------------
% Read input arguments
if ~iscell(mean)
    MU = mean;
else
    if numel(mean) >= 1
        MU = mean{1};
        if numel(mean) >= 2
            b = mean{2};
        end
    end
end
if ~iscell(prec)
    V = prec;
else
    if numel(prec) >= 1
        V = prec{1};
        if numel(prec) >= 2
            n = prec{2};
        end
    end
end


%--------------------------------------------------------------------------
% Use double precision
MU = double(MU);
b  = double(b);
V  = double(V);
n  = double(n);

%--------------------------------------------------------------------------
% Dimensions
P = size(MU,1);
K = size(MU,2);
    
if nargin == 3 && ~isempty(L)
%--------------------------------------------------------------------------
% Marginal distribution (do not use missing dimensions)

% Assume that:
% . A is a KxK positive-definite matrix
% . A ~ W_K(V,n)
% . S = inv(A)
% . A is partitioned as [A11 A12: A12' A22] with A11 PxP
% . V is partitioned as [V11 V12: V12' V22] with V11 PxP
% . S is partitioned as [S11 S12: S12' S22] with S11 PxP
% Then
% . A11 ~ W_P(V11,n)
% . inv(S11) = A11 - A12*A22\A12'
% . inv(S11) ~ W_P(V11 - V12*V22\V12', n - K + P)
% This allows us to compute E[inv(S11)] and E[ln|S11]], which are needed to
% compute the expected marginal distribution within each cluster.

    c = zeros(numel(L), K, 'like', MU);
    for i=1:numel(L)
        code = L(i);
        observed = code2bin(code, P);
        missing  = ~observed;
        Pm = sum(missing);
        Po = sum(observed);
        for k=1:K
            Vo = V(observed,observed,k) - V(observed,missing,k)*(posdef.inv(V(missing,missing,k))*V(missing,observed,k));
            c(i,k) = - 0.5 * Po * log(2*pi);
            if sum(n) > 0
                no = n(k) - Pm;
                c(i,k) = c(i,k) + 0.5 * wishart.Elogdet(Vo,no) ...
                                - 0.5 * no * MU(observed,k).' * Vo * MU(observed,k);
            else
                c(i,k) = c(i,k) + 0.5 * posdef.logdet(Vo) ...
                                - 0.5 * MU(observed,k).' * Vo * MU(observed,k);
            end
            if sum(b) > 0
                c(i,k) = c(i,k) - 0.5 * Po / b(k);
            end
        end
        
    end
    
else
%--------------------------------------------------------------------------
% No missing dimensions
    c = zeros(1,K, 'like', MU);
    for k=1:K
        c(k) = - 0.5 * P * log(2*pi);
        if sum(n) > 0
            c(k) = c(k) + 0.5 * wishart.Elogdet(V(:,:,k),n(k)) ...
                        - 0.5 * n(k) * MU(:,k)' * V(:,:,k) * MU(:,k);
        else
            c(k) = c(k) + 0.5 * posdef.logdet(V(:,:,k)) ...
                        - 0.5 * MU(:,k)' * V(:,:,k) * MU(:,k);
        end
        if sum(b) > 0
            c(k) = c(k) - 0.5 * P / b(k);
        end
    end
end