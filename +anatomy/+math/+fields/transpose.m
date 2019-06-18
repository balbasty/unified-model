function F = transpose(F, dims)
%__________________________________________________________________________
% anatomy.math.fields.transpose
%--------------------------------------------------------------------------
% FORMAT F = transpose(F, [dims])
% F    - A field of matrices (or vectors)
% dims - Dimensions of the nd-array corresponding to the matrix dimensions
%           default: [2 3]
%
% Non-conjugate transpose of a field of matrices.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

    % By default: assume a three-dimensional field
    if nargin < 2 || isempty(dims)
        dims = [2 3];
    end
    % If only one dimension provided, assum we permute with the next one
    if numel(dims) == 1
        dims = [dims dims+1];
    end
    % If more than two dimensions provided, issue a warning
    if numel(dims) > 2
        warning('Only two dimensions required. I will discard the remaining ones')
        dims = dims(1:2);
    end

    perms = 1:max(numel(size(F)), max(dims));
    perms(dims(1)) = dims(2);
    perms(dims(2)) = dims(1);
    F = permute(F, perms);

end