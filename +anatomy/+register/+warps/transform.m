function y = transform(A, y)
%__________________________________________________________________________
% anatomy.register.warps.transform
%--------------------------------------------------------------------------
% FORMAT y = transform(A, y)
% iy - A warp (Nx * Ny * Nz * 3) 
% A  - an affine transformation matrix (4 * 4)
% oy - A warp with the same lattice a iy.
%
% This is equivalent to composing the affine transform and the warp: A o iy
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging
    dim = size(y);
    y   = reshape(y, [prod(dim(1:3)), 3]);
    y   = bsxfun(@plus, y*A(1:3,1:3)', A(1:3,4)');
    y   = reshape(y, dim);
end