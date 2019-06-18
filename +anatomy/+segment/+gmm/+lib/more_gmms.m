function gmm = more_gmms(gmm,part)
%__________________________________________________________________________
% anatomy.segment.gmm.lib.more_gmms
%--------------------------------------------------------------------------
% FORMAT gmm = gmm.lib.more_gmms(gmm, part)
%
% gmm  - Cell with the following format {m,b,W,n}, where there are K
%        Gaussians.
% part - [1,K_p] vector partitioning a K GMM into a K_p GMM (K_P>=K). E.g.
%        [1 1 1 2 3 4 5 6 6] means that the first Gaussian will be divided
%        into 3 and the last into 2. The rest will remain the same.
%
% gmm  - Cell with the following format {m,b,W,n}, where there are K_p 
%        Gaussians.
%
% A crude heuristic to replace a single VB Gaussian by a bunch of VB Gaussians.
% If there is only one Gaussian, then it should be the same as the
% original distribution.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

MU0 = gmm{1};
b0  = gmm{2};
W0  = gmm{3};
n0  = gmm{4};

K   = numel(n0);
K_p = numel(part);
C   = size(MU0,1);

A0 = bsxfun(@times,W0,reshape(n0,[1 1 K])); % E[Lambda]

m = zeros(C,K_p);
b = zeros(1,K_p);
W = zeros(C,C,K_p);
n = zeros(1,K_p); % A = nW, W = 1/n*inv(Cov)

for k=1:K    
    kk = sum(part==k);
    w  = 1./(1 + exp(-(kk - 1)*0.25)) - 0.5;
    mn = MU0(:,k);
    vr = inv(A0(:,:,k));
    
    mn = sqrtm(vr)*randn(C,kk)*w + repmat(mn,[1 kk]);
    vr = vr*(1 - w);
    pr = inv(vr);
    W1 = (1/n0(k))*pr;
    
    m(:,part==k)   = mn;
    b(part==k)     = b0(k);
    W(:,:,part==k) = repmat(W1,[1 1 kk]);
    n(part==k)     = n0(k);
end

gmm{1} = m;
gmm{2} = b;
gmm{3} = W;
gmm{4} = n;