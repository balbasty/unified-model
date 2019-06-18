function [klMU,klA] = gausswishart(mean,prec,mean0,prec0)
%__________________________________________________________________________
% anatomy.segment.gmm.lib.kl.gausswishart
%--------------------------------------------------------------------------
% [klMU,klA] = gmm.lib.kl.gausswishart({MU,b}, {V,n}, {MU0,b0}, {V0,n0})
%
% KL divergence between two Gauss-Wishart distributions
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

import anatomy.math.matrix.*       % posdef.logdet
import anatomy.math.probability.*  % wishart.Elogdet, wishart.kl

MU  = [];
b   = [];
V   = []; % It can actually be A (when n == 0)
n   = [];
MU0 = [];
b0  = [];
V0  = [];
n0  = [];

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
if nargin >= 3
    if ~iscell(mean0)
        MU0 = mean0;
    else
        if numel(mean0) >= 1
            MU0 = mean0{1};
            if numel(mean0) >= 2
                b0 = mean0{2};
            end
        end
    end
end
if nargin >=4
    if ~iscell(prec0)
        V0 = prec0;
    else
        if numel(prec0) >= 1
            V0 = prec0{1};
            if numel(prec0) >= 2
                n0 = prec0{2};
            end
        end
    end
end

%--------------------------------------------------------------------------
% Read input arguments
P = size(MU,1);
K = size(MU,2);
LogDetA = zeros(1,K, 'like', V);
if sum(n) > 0
    A = bsxfun(@times, V, reshape(n, [1 1 K]));
    for k=1:K
        LogDetA(k) = wishart.Elogdet(V(:,:,k),n(k));
    end
else
    A = V;
    for k=1:K
        LogDetA(k) = posdef.logdet(A(:,:,k));
    end
end

% Lower bound
klMU = zeros(1, 'like', MU);
klA  = zeros(1, 'like', A);
for k=1:K
    % + prior
    if sum(b0) > 0
        % prior
        klMU = klMU - P*log(2*pi) ...
                    + P*log(b0(k)) ...
                    + LogDetA(k) ...
                    - b0(k)*(MU(:,k)-MU0(:,k)).'*A(:,:,k)*(MU(:,k)-MU0(:,k)) ...
                    - P*b0(k)/b(k);
        % posterior
        klMU = klMU + P*log(2*pi) ...
                    - P*log(b(k)) ...
                    - LogDetA(k) ...
                    + P;
    end
    if sum(n0) > 0
        klA = klA - wishart.kl(V(:,:,k), n(k), V0(:,:,k), n0(k));
    end
end
klMU = 0.5 * klMU;