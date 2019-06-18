function [Lam1, n1] = update(varargin)
%__________________________________________________________________________
% anatomy.math.probability.wishart.update
%--------------------------------------------------------------------------
% FORMAT [Lam1, n1]  = update(Lam, n,        Lam0, n0)
% FORMAT [Lam1, n1]  = update(ss0, ss1, ss2, Lam0, n0, (mu=0))
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging
    
    import anatomy.math.matrix                 % posdef.inv
    
    if nargin > 4
        % Sufficient statistics case
        ss0  = varargin{1};
        ss1  = varargin{2};
        ss2  = varargin{3};
        Lam0 = varargin{4};
        n0   = varargin{5};
        if nargin < 6
            iV = ss2;
        else
            iV = ss2 - 2*ss1*mu' + ss0*mu*mu';
        end
        n1   = n0+ss0;
        if n0, Lam1 = (n0*posdef.inv(Lam0) + iV)/n1;
        else,  Lam1 = iV/n1; end
        % stable inverse
        Lam1 = posdef.inv(Lam1);
    else
        % Average case
        Lam  = varargin{1};
        n    = varargin{2};
        Lam0 = varargin{3};
        n0   = varargin{4};
        n1   = n + n0;
        if n0, Lam1 = n1*posdef.inv(...
                            n0 * posdef.inv(Lam0) + ...
                            n  * posdef.inv(Lam));
        else,  Lam1 = Lam; end
    end
    
    
end