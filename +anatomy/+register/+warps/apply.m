function yf = apply(f, y, itrp, bnd)
%__________________________________________________________________________
% anatomy.register.warps.apply
%--------------------------------------------------------------------------
% FORMAT yf = apply(f, y, (itrp), (bnd))
% f     - Input image (Nx * Ny * Nz * ...)
% y     - Non-linear warp (Mx * My * Mz * 3)
% itrp  - Interpolation order [default: 1 1 1]
% bnd   - Boundary conditions (0/1 = mirror/circulant) [default: 1 1 1]
% yf    - Deformed image  (Mx * My * Mz * ...)
%
% Warps an image with a non-linear transform, i.e., computes f(y).
% The input image can be non-scalar.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

    if nargin < 4
        bnd = [1 1 1];
        if nargin < 3
            itrp = [1 1 1];
        end
    end
    
    if numel(itrp) < 3
        itrp = padarray(itrp, [0 3 - numel(itrp)], 'replicate', 'post');
    end
    itrp = double(itrp(1:3));
    itrp = itrp(:)';
    if numel(bnd) < 3
        bnd = padarray(bnd, [0 3 - numel(bnd)], 'replicate', 'post');
    end
    bnd = double(bnd(1:3));
    bnd = bnd(:)';
    
    dim_in = size(f);
    f = reshape(f, dim_in(1), dim_in(2), dim_in(3), []);
    dim_out = size(y);
    dim_out = dim_out(1:3);
    
    yf = zeros([dim_out size(f, 4)], 'like', f);
    for k=1:size(f, 4)
        yf(:,:,:,k) = apply_scalar(f(:,:,:,k), y, itrp, bnd);
    end
    yf = reshape(yf, [dim_out dim_in(4:end)]);
    
    
function yf = apply_scalar(f, y, itrp, bnd)
% FORMAT out = warp_scalar(in, y,  (itrp))
% in    - Input _scalar_ image (or function R^3 -> R).
% y     - Non-linear warp.
% itrp  - Interpolation order.
% bnd   - Boundary conditions (0/1 = mirror/circulant).
%
% Warps an image with a non-linear transform, i.e., computes in(y).
% The input image must be scalar.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

    % Interpolate input image with circulant boundaries
    coeff = spm_diffeo('bsplinc', single(f), [itrp bnd]);
    
    % Interpolate image on output grid
    yf = spm_diffeo('bsplins', coeff, single(y), [itrp bnd]);