function code = obs2code(X)
%__________________________________________________________________________
% anatomy.segment.gmm.lib.obs2code
%--------------------------------------------------------------------------
% FORMAT code = gmm.lib.obs2code(X)
%
% Compute a "missing code" image for the input observation matrix.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

import anatomy.segment.gmm.lib.double2int

code = double2int(sum(bsxfun(@times, ~isnan(X), 2.^(0:size(X,2)-1)), 2));