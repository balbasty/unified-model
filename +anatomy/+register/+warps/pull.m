function yf = pull(f, y, bnd, xtr)
%__________________________________________________________________________
% anatomy.register.warps.pull
%--------------------------------------------------------------------------
% FORMAT yf = pull(f, y, (bnd), (xtr))
% f     - Input image (Nx * Ny * Nz * ...)
% y     - Non-linear warp (Mx * My * Mz * 3)
% bnd   - Boundary conditions (0/1 = circulant/neumann) [default: 0]
% xtr   - Extrapolate data outside of the FOV [false]
% yf    - Pulled image  (Mx * My * Mz * ...)
%
% Pull an image with a non-linear transform, i.e., computes f(y).
% This is closely related to warps.apply. However, here the interpolation
% order is restricted to 1, *but* data can be extrapolated with any
% boundary conditions.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

    if nargin < 4
        xtr = false;
        if nargin < 3
            bnd = 0;
        end
    end
    
    if ischar(bnd)
        bnd = double(lower(bnd(1)) == 'n');
    end
    bnd = double(bnd(1));
    xtr = double(xtr(1));
    
    dim_in = size(f);
    f = reshape(f, dim_in(1), dim_in(2), dim_in(3), []);
    dim_out = size(y);
    dim_out = dim_out(1:3);
    
    yf = zeros([dim_out size(f, 4)], 'like', f);
    for k=1:size(f, 4)
        yf(:,:,:,k) = pull_scalar(f(:,:,:,k), y, bnd, xtr);
    end
    yf = reshape(yf, [dim_out dim_in(4:end)]);
    
    
function yf = pull_scalar(f, y, bnd, xtr)
    
    if xtr
        func = 'pullc';
    else
        func = 'pull';
    end

    % Set boundary conditions
    spm_diffeo('boundary', bnd);
    % Interpolate image on output grid
    yf = spm_diffeo(func, single(f), single(y));