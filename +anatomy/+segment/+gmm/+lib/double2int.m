function L = double2int(L)
%__________________________________________________________________________
% anatomy.segment.gmm.lib.double2int
%--------------------------------------------------------------------------
% FORMAT L = gmm.lib.double2int(L)
%
% Find the best suited integer type to convert L, based on min and max
% values
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

import anatomy.segment.gmm.lib.range2int

minval = min(L(:));
maxval = max(L(:));
type   = range2int(maxval,minval);
func   = str2func(type);
L      = func(L);