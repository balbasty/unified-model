function logpX = marginal(X, cluster, const, codes, E)
%__________________________________________________________________________
% anatomy.segment.gmm.lib.marginal
%--------------------------------------------------------------------------
% logp = gmm.lib.marginal(X, {MU,A}, const, {C, L}, E)
% logp = gmm.lib.marginal(X, {MU,V,n}, const, {C, L}, E)
% 
% X         - NxP   Observed values
% MU        - PxK   (Expected) means
% A         - PxPxK (Expected) precision matrices
% const     - MxK   Constant terms. If M > 1, marginal distributions.
% C         - Nx1   Image of "missing" codes
% L         - Mx1   List of existing codes
% E         - 1xP   Binning uncertainty
%
% logpX     - NxK   (Expected) log-likelihood of belonging to each class
%
% Compute the expected log-likelihood of each observation belonging to each
% cluster: logpx(i,k) = E[ln p(g(i) | MU_k,A_k)]
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

import anatomy.segment.gmm.lib.code2bin
import anatomy.math.matrix.*               % posdef.inv

MU = [];
A  = [];
V  = [];
n  = [];
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
            if numel(cluster) >= 3
                V = A;
                n = cluster{3};
            end
        end
    end
end
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
if nargin < 5
    E = [];
end

% -------------------------------------------------------------------------
% Use double precision
X  = double(X);
MU = double(MU);
A  = double(A);
V  = double(V);
n  = double(n);
E  = double(E);
const = double(const);

% -------------------------------------------------------------------------
% Dimensions
N  = size(X,1);
P  = size(X,2);
K  = size(MU,2);
if isempty(L), L = 2^P - 1;    end % None missing
if isempty(E), E = zeros(1,P, 'like', X); end % No uncertainty
if size(const, 1) == 1
    const = repmat(const, [numel(L) 1]);
end

logpX  = zeros([N K], 'like', X);

% -------------------------------------------------------------------------
% For each combination of missing voxels
for i=1:numel(L)
    
    % ---------------------------------------------------------------------
    % Get mask of missing values and modalities (with this particular code)
    c        = L(i);
    observed = code2bin(c, P);
    Po       = sum(observed);
    if Po == 0, continue; end
    if isempty(C) && (c == 2^P-1),  msk = ones(N,1, 'logical');
    else,                           msk = (C == c);   end
    missing  = ~observed;
    Pm       = P-Po;
    Nm       = sum(msk(:));
    if Nm == 0, continue; end
    
    % ---------------------------------------------------------------------
    % Initialise with constant term
    logpX(msk,:) = repmat(const(i,:), [Nm 1]);
    
    % ---------------------------------------------------------------------
    % Non constant terms
    X1 = X(msk,observed)';
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
        
        % Quadratic term in observed values: (obs-mean) x (obs-mean)
        l = Ao * bsxfun(@minus, X1, 2*MU(observed,k));
        l = -0.5 * dot(l, X1, 1);
        
        % Binning uncertainty
        if any(any(E))
            if size(E,1)==1
                l = l - 0.5 * trace(Ao * diag(E(observed)));
            else
                l = l - 0.5 * sum(bsxfun(@times,diag(Ao)',E(msk,observed)),2)';
            end
        end
        
        % Reshape as a column vector
        logpX(msk,k) = logpX(msk,k) + l';
    end
        
end

% Set NaN value for voxels without observed dimensions
logpX(all(isnan(X),2),:) = NaN;