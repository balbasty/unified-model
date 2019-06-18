function F = ctranspose(F, varargin)
%__________________________________________________________________________
% anatomy.math.fields.ctranspose
%--------------------------------------------------------------------------
% FORMAT F = ctranspose(F, [dims])
% F    - A field of matrices (or vectors)
% dims - Dimensions of the nd-array corresponding to the matrix dimensions
%           default: [2 3]
%
% Conjugate transpose of a field of matrices.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

    import anatomy.math.fields.transpose
    F = conj(transpose(F, varargin{:}));
    
end