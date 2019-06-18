function id = identity(dim)
%__________________________________________________________________________
% anatomy.register.warps.identity
%--------------------------------------------------------------------------
% id = identity(dim)
%
% Generate an identity warp for a lattice of dimensions `dim`.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging
    id        = cell(3,1);
    [id{1:3}] = ndgrid(single(1:dim(1)),single(1:dim(2)),single(1:dim(3)));
    id        = cat(4,id{:});
end