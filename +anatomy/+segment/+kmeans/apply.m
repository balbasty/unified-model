function [L,D] = apply(X, C, dist, missing)
%__________________________________________________________________________
% anatomy.segment.kmeans.apply
%--------------------------------------------------------------------------
% FORMAT [L,D] = kmeans.apply(X, C)
% Classify observations based on known centroids.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

dist = str2func(dist);

% Reshape input
dim = size(X);
if dim(end) == size(C,2)
    dim = dim(1:end-1);
end
X = reshape(X, [], size(C,2));

% Compute distance
D = dist(X,C);
% Get closest centroid
[~, L] = min(D, [], 2);
L      = single2int(L);
% Replace missing data
if ~missing
    msk      = any(isnan(X),2);
    L(msk)   = 0;
    D(msk,:) = NaN; 
end
% Reshape output
L = reshape(L, [dim 1]);
D = reshape(D, [dim size(C,1)]);