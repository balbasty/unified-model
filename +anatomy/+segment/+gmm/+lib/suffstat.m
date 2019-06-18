function varargout = suffstat(varargin)
%__________________________________________________________________________
% anatomy.segment.gmm.lib.suffstat
%--------------------------------------------------------------------------
% FORMAT [SS0,SS1,SS2] = gmm.lib.suffstat(X, Z, W)
% FORMAT [SS0,SS1,SS2] = gmm.lib.suffstat('base',  X, Z, W, {C,L})
% FORMAT [SS0,SS1,SS2] = gmm.lib.suffstat('infer', SS0, SS1, SS2, {MU,A}, L)
% FORMAT         [SS2] = gmm.lib.suffstat('bin',   E, Z, W, {C,L})
%
% X    - NxP Observed + Inferred values
% S    - Posterior uncertainty for each code Cx{Nmx(Km(Km+1)/2)}
% W    - Nx1 Observation weights
% E    - 1xP Binning uncertainty
% Z    - NxK Responsibilities
% C    - Nx1 Missing values "code"
% L    - List of codes present (saves a tiny bit of time if provided)
%
% SS0 - 1xK   0th order suff stat (sum of resp)
% SS1 - PxK   1st order suff stat (weighted sum of intensities)
% SS2 - PxPxK 2nd order suff stat (weighted sum of squared intensities)
%
% Compute sufficient statistics up to 2nd order, taking into account
% inferred values and their uncertainty.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

if nargin == 0
    help suffstat
    error('Not enough argument. Type ''suffstat'' for help.');
end
if ~ischar(varargin{1})
    [varargout{1:nargout}] = suffstat_default(varargin{:});
    return
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case {'base'}
        [varargout{1:nargout}] = suffstat_base(varargin{:});
    case {'infer'}
        [varargout{1:nargout}] = suffstat_infer(varargin{:});
    case {'bin'}
        [varargout{1:nargout}] = suffstat_bin(varargin{:});              
    otherwise
        help suffstat
        error('Unknown function %s. Type ''help suffstat'' for help.', id)
end

% =========================================================================
function [SS0,SS1,SS2] = suffstat_default(X, Z, W)
% FORMAT [SS0,SS1,SS2] = suffstat_default(X, Z, W)
%
% Compute sufficient statistics (up to 2nd order)

%--------------------------------------------------------------------------
% Dimensions
N = size(X, 1);
P = size(X, 2);
K = size(Z, 2);

%--------------------------------------------------------------------------
% Use double precision
X = double(X);
Z = double(Z);
W = double(W);

%--------------------------------------------------------------------------
% Weight responsibilities
Z = bsxfun(@times, Z, W);

%--------------------------------------------------------------------------
% Oth order
SS0 = sum(Z, 1, 'omitnan');

%--------------------------------------------------------------------------
% 1st order
SS1 = sum(bsxfun(@times, X, reshape(Z, [N 1 K])), 1, 'omitnan');
SS1 = reshape(SS1, [P K]);

%--------------------------------------------------------------------------
% 1nd order
SS2 = zeros(P,P,K, 'like', Z);
for i=1:P
    SS2(i,i,:) = reshape(sum(bsxfun(@times, Z, X(:,i).^2),1,'omitnan'), [1 1 K]);
    for j=i+1:P
        SS2(i,j,:) = reshape(sum(bsxfun(@times, Z, X(:,i).*X(:,j)),1,'omitnan'), [1 1 K]);
        SS2(j,i,:) = SS2(i,j,:);
    end
end

% =========================================================================
function [SS0,SS1,SS2] = suffstat_base(X, Z, W, codes)
% FORMAT [{SS0},{SS1},{SS2}] = suffstat_base(X, Z, W, {C,L})
%
% Compute sufficient statistics (up to 2nd order)

import anatomy.segment.gmm.lib.code2bin

if nargin < 3
    W = 1;
end


C  = [];
L  = [];

%--------------------------------------------------------------------------
% Read input arguments
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
% Use double precision
X = double(X);
Z = double(Z);
W = double(W);

%--------------------------------------------------------------------------
% Dimensions
N = size(X,1);
P = size(X,2);
K = size(Z,2);
Z = bsxfun(@times, Z, W); % Multiply resp with observation count
if isempty(L)
    L = 2^P-1;
end

SS0 = cell(1,numel(L));
SS1 = cell(1,numel(L));
SS2 = cell(1,numel(L));

%--------------------------------------------------------------------------
% Sum missing data
for i=1:numel(L)
    c       = L(i);
    missing = ~code2bin(c, P);
    if isempty(C),  msk = ones(N, 1, 'logical');
    else,           msk     = (C == c);
    end
    Pm      = sum(missing);
    Po      = P-Pm;
    Nm      = sum(msk);
    if Nm == 0, continue; end

    X1 = X(msk,~missing);
    Z1 = Z(msk,:);
    
    % ---------------------------------------------------------------------
    % Oth order moment
    SS0{i} = sum(Z1, 1);
    if nargout == 1, return; end

    % ---------------------------------------------------------------------
    % 1st order moment
    SS1{i} = sum(bsxfun(@times, X1, reshape(Z1, [Nm 1 K])), 1);
    SS1{i} = reshape(SS1{i}, [Po K]);
    if nargout == 2, return; end

    % ---------------------------------------------------------------------
    % 2nd order moment
    SS2{i} = zeros(Po,Po,K, 'like', Z);
    for k=1:Po
        SS2{i}(k,k,:) = reshape(sum(bsxfun(@times, Z1, X1(:,k).^2),1), [1 1 K]);
        for j=k+1:Po
            SS2{i}(k,j,:) = reshape(sum(bsxfun(@times, Z1, X1(:,k).*X1(:,j)),1), [1 1 K]);
            SS2{i}(j,k,:) = SS2{i}(k,j,:);
        end
    end
end


% =========================================================================
function [SS0,SS1,SS2] = suffstat_infer(lSS0, lSS1, lSS2, cluster, L)
% FORMAT [SS0,SS1,SS2] = suffstat_infer(SS0, SS1, SS2, {MU,A}, L)
%
% lSS0 - List of sufficient statistics for each pattern of missing data
% lSS1 - List of sufficient statistics for each pattern of missing data
% lSS2 - List of sufficient statistics for each pattern of missing data
% MU   - Clusters' mean
% A    - Clusters' precision matrix
% L    - List of unique codes
% 
% Compute "missing" 1st/2nd order statistics based on a simplified
% inferrence of missing values.
% simplified = inference is performed cluster-wise

import anatomy.segment.gmm.lib.code2bin
import anatomy.math.matrix.*               % posdef.inv

MU = [];
A  = [];

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
if nargin < 5
    L = [];
end

%--------------------------------------------------------------------------
% Use double precision
MU = double(MU);
A  = double(A);

%--------------------------------------------------------------------------
% Dimensions
P = size(MU,1);
K = size(MU,2);

SS0 = zeros(1,K, 'like', lSS0{1});
if nargout > 1
    SS1 = zeros(P,K, 'like', lSS1{1});
    if nargout > 2
        SS2 = zeros(P,P,K, 'like', lSS2{1});
    end
end

for i=1:numel(L)
    c        = L(i);
    observed = code2bin(c, P);
    missing  = ~observed;
    
    for k=1:K
        % -----------------------------------------------------------------
        % 0th order moment
        SS0k   = lSS0{i}(k);
        SS0(k) = SS0(k) + SS0k;
        
        if nargout > 1
        % -----------------------------------------------------------------
        % 1st order moment
            SS1k = lSS1{i}(:,k);
            MUo  = MU(observed,k);
            MUm  = MU(missing,k);
            SA   = posdef.inv(A(missing,missing,k)) * A(missing,observed,k);
            
            % 1) observed
            SS1(observed,k) = SS1(observed,k) + SS1k;
        
            % 2) missing
            % > t = mu(m) + A(m,m) \ A(m,o) * (mu(o) - g)
            SS1(missing,k) = SS1(missing,k) + SS0k * MUm;
            SS1(missing,k) = SS1(missing,k) + SA * (SS0k * MUo - SS1k);
        end
        
        if nargout > 2
        % -----------------------------------------------------------------
        % 2nd order moment: quadratic terms
            SS2k = lSS2{i}(:,:,k);
        
            % 0) precompute stuff
            MUMUm = SS0k * (MUm * MUm.');
            MUMUo = SS0k * (MUo * MUo.');
            GMUo  = SS1k * MUo.';
            GMUm  = SS1k * MUm.';
        
            % 1) observed x observed
            SS2(observed,observed,k) = SS2(observed,observed,k) + SS2k;
            
            % 2) missing x observed
            tmp = GMUm.' + SA * (GMUo.' - SS2k);
            SS2(missing,observed,k) = SS2(missing,observed,k) + tmp;
            SS2(observed,missing,k) = SS2(observed,missing,k) + tmp.';
            
            % 3) missing x missing
            SS2(missing,missing,k) = SS2(missing,missing,k) + MUMUm;
            tmp = SA * (SS0k * MUo - SS1k) * MUm.';
            SS2(missing,missing,k) = SS2(missing,missing,k) + tmp + tmp';
            SS2(missing,missing,k) = SS2(missing,missing,k) ...
                + SA * (SS2k + MUMUo - GMUo.' - GMUo) * SA.';
    
            % 4) uncertainty ~ missing
            SS2(missing,missing,k) = SS2(missing,missing,k) ...
                + SS0k * podef.inv(A(missing,missing,k));
        end
    end
end


% =========================================================================
function SS2 = suffstat_bin(E, Z, W, codes)
% FORMAT SS2 = suffstat_bin(E, Z, W, {C,L})
%
% E - Variance in each modality due to binning
% Z - Responisbilities
% W - Observation weights
% C - Missing code image
% L - List of unique codes
% P - Observed space dimension
% 
% Compute "uncertainty" 2nd order statistics based on the posterior
% precision matrix about inferred values.

import anatomy.segment.gmm.lib.code2bin

C  = [];
L  = [];
if nargin < 4
    if nargin < 3
        W = 1;
    end
end

%--------------------------------------------------------------------------
% Read input arguments
if nargin >= 5
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
% Use double precision
E = double(E);
Z = double(Z);
W = double(W);

%--------------------------------------------------------------------------
% Dimensions
N = size(Z,1);
K = size(Z,2);
P = size(E,2);
if sum(E) == 0
    SS2 = zeros(P,P,size(Z,2), 'like', Z);
    return
end
if isempty(L)
    L = 2^P - 1; % None missing
end
Z   = bsxfun(@times, Z, W); % Multiply resp with observation count
SS2 = zeros(P,P,K, 'like', Z);

% -------------------------------------------------------------------------
% 2nd order moment: uncertainty ~ binning
for i=1:numel(L)
    c        = L(i);
    observed = code2bin(c, P);
    Pp       = sum(observed);
    if Pp == 0, continue; end
    
    if isempty(C) && (c == 2^P-1),  msk = ones(N,1,'logical');
    else,                           msk = (C == c);  end
    Nm = sum(msk);
    if Nm == 0, continue; end
    
    
    list_p = 1:P;
    list_p = list_p(observed);
    for p=list_p
        if size(E,1)==1
            SS2(p,p,:) = SS2(p,p,:) ...
                + bsxfun(@times, reshape(sum(Z(msk,:), 1), [1 1 K]), E(p));
        else
            SS2(p,p,:) = SS2(p,p,:) ...
                + reshape(sum(bsxfun(@times, Z(msk,:), E(msk,p)), 1), [1 1 K]);
        end
    end
end