function [lb,cst] = marginalsum(SS0, SS1, SS2, mean, prec, L, SS2b)
%__________________________________________________________________________
% anatomy.segment.gmm.lib.marginalsum
%--------------------------------------------------------------------------
% [lb,const] = gmm.lib.marginalsum(SS0, SS1, SS2, MU, A, L, SS2b)
% [lb,const] = gmm.lib.marginalsum(SS0, SS1, SS2, {MU,b}, {V,n}, L, SS2b)
% 
% SS0       - {1xK}   Zero-th order moment (per config)
% SS1       - {PxK}   First   order moment (per config)
% SS1       - {PxPxK} Second  order moment (per config)
% MU        - PxK     Means
% A/V       - PxPxK   Precision/Scale matrices
% b         - 1xK     Mean degrees of freedom
% n         - 1xK     Precision degrees of freedom
% L         - Mx1     List of existing codes
% SS2b      - PxPxK   Binning uncertainty
%
% lb        -         Sum of (expected) marginal likelihoods
% const     - MxK     Constant terms
%
% Compute the expected log-likelihood of each observation belonging to each
% cluster: lb = sum_{i,k} E[z_ik] E[ln p(g(i) | MU_k,A_k)]
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

import anatomy.segment.gmm.lib.code2bin
import anatomy.math.matrix.*               % posdef.inv

MU = [];
A  = [];
V  = [];
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
    A = prec;
else
    if numel(prec) >= 1
        A = prec{1};
        if numel(prec) >= 2
            V = A;
            n = prec{2};
        end
    end
end
if nargin <= 6
    L = [];
end
if nargin <= 7
    SS2b = 0;
end

% -------------------------------------------------------------------------
% Dimensions
P  = size(MU,1);
K  = size(MU,2);
if isempty(L), L = 2^P - 1;    end % None missing

% -------------------------------------------------------------------------
% Constant term
cst = const(mean, prec, L);

lb  = 0;

% -------------------------------------------------------------------------
% For each combination of missing voxels
for i=1:numel(L)
    
    % ---------------------------------------------------------------------
    % Get mask of missing values and modalities (with this particular code)
    c        = L(i);
    observed = code2bin(c, P);
    Po       = sum(observed);
    if Po == 0, continue; end
    missing  = ~observed;
    Pm       = P-Po;
    
    % ---------------------------------------------------------------------
    % Initialise with constant term
    lb = lb + sum(cst(i,:) .* SS0{i});
    
    % ---------------------------------------------------------------------
    % Non constant terms
    for k=1:K
        
        % /!\ Sub-covariance is different from the inverse sub-precision
        % ML case:
        %   inv(S(o,o)) = A(o,o) - A(o,m)*A(m,m)\A(m,o)
        % Bayesian case:
        %   inv(S(o,o)) ~ W(V(o,o) - V(o,m)*V(m,m)\V(m,o), n - Pm)
        if sum(n) > 0
            Ao = V(observed,observed,k) - V(observed,missing,k)*(posdef.inv(V(missing,missing,k))*V(missing,observed,k));
            Ao = (n(k)-Pm) * Ao;
        else
            Ao = A(observed,observed,k) - A(observed,missing,k)*(posdef.inv(A(missing,missing,k))*A(missing,observed,k));
        end
        
        % 1) obs x mean
        lb = lb + SS1{i}(:,k).' * Ao * MU(observed,k);
        
        % 1) obs x obs
        lb = lb - 0.5 * trace(Ao * SS2{i}(:,:,k));
        
        % 3) Binning uncertainty
        lb = lb - 0.5 * trace(Ao * SS2b);
        
    end
        
end