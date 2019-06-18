function A = inv(A)
%__________________________________________________________________________
% anatomy.math.matrix.posdef.inv
%--------------------------------------------------------------------------
% FORMAT iA = inv(A)
% A  - A positive-definite square matrix
% iA - Its inverse
%
% Stable inverse of a positive-definite matrix.
% Eigendecomposition is used to compute a more stable inverse.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

% John Ashburner

    % Import loaddiag function
    import anatomy.math.matrix.loaddiag

    % Use eigendecomposition to invert
    [V,D] = eig(A);
    if any(diag(D) <= 0)
        warning('Matrix has negative eigenvalues')
        D(D <= 0) = eps; % Threshold negative eigenvalues
    end
    D     = loaddiag(D);
    A     = real(V * (D \ V'));

end