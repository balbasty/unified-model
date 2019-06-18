function A = loaddiag(A)
%__________________________________________________________________________
% anatomy.math.matrix.loaddiag
%--------------------------------------------------------------------------
% FORMAT A = loaddiag(A)
% A  - A square matrix
%
% Load A's diagonal until it is well conditioned for inversion.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

    factor = 1e-7;
    while rcond(A) < 1e-5
        A = A + factor * max([diag(A); eps]) * eye(size(A));
        factor = 10 * factor;
    end

end