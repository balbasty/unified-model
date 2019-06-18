function y = compose(varargin)
%__________________________________________________________________________
% anatomy.register.warps.compose
%--------------------------------------------------------------------------
% FORMAT y = compose(y_1, ..., y_n, (itrp))
% y_i  - A warp OR an affine matrix.
% itrp - Interpolation degree [default: 1]
% y    - A warp with the same lattice and voxel size as y_n.
%
% Compose a series of transformations.
% NB:
% - The right-most transform should always be a warp.
% - To specify the output lattice, add an identity warp as the
%   right-most argument.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging
    
    import anatomy.register.warps.transform
    import anatomy.register.warps.identity
    import anatomy.register.warps.apply

    % Get interpolation order
    if ~isempty(varargin) && isscalar(varargin{end})
        itrp     = varargin{end};
        varargin = varargin(1:end-1);
    else
        itrp = 1;
    end
    
    if isempty(varargin)
        y = [];
        return
    end
    
    % Loop
    while ~isempty(varargin)
        cur_y = varargin{end};
        varargin = varargin(1:end-1);
        if length(size(cur_y)) == 2
            % Affine o Warp
            y = transform(cur_y, y);
        else
            % Warp o Warp
            % > Warp the deformation field with circulant boundary
            cur_lat = size(cur_y);
            cur_lat = cur_lat(1:3);
            id      = identity(cur_lat);
            y       = apply(cur_y - id, y, itrp, 1) + y;
        end
    end
end