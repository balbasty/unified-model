function M = eigthreshold(M)
%__________________________________________________________________________
% anatomy.math.matrix.eigthreshold
%--------------------------------------------------------------------------
% FORMAT M = eigthreshold(M)
%
% Stabilise matrix by thresholding negative eigenvalues
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

    [V,D1] = eig(M);
    tol    = max(diag(D1))*eps('single');
    D1     = diag(max(diag(D1),tol));
    M      = real(V*D1*V');
end