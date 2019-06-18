function el = Elog(varargin)
%__________________________________________________________________________
% anatomy.math.probability.wishart.beta.Elog
%--------------------------------------------------------------------------
% FORMAT el = Elog(a, b)
% FORMAT el = Elog(p, n, 'ber')
% FORMAT el = Elog(p, n, k, 'bin')
% FORMAT el = Elog(p, n, r, 'nbin')
% FORMAT el = Elog(p, n, 'geom')
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    import anatomy.math.probability.beta   % Elog
    
    % Check if we are in the reparameterised case
    if nargin > 2
        switch lower(varargin{end})
            case 'ber'
                p = varargin{1};
                n = varargin{2};
                a = bsxfun(@times, n, p);
                b = bsxfun(@times, n, 1-p);
            case 'bin'
                p = varargin{1};
                n = varargin{2};
                k = varargin{3};
                a = bsxfun(@times, k, bsxfun(@times, n, p));
                b = bsxfun(@times, k, bsxfun(@times, n, 1-p));
            case 'nbin'
                p = varargin{1};
                n = varargin{2};
                r = varargin{3};
                b = bsxfun(@times, r, n);
                a = bsxfun(@rdivide, bsxfun(@times, p, b), 1-p);
            case 'geom'
                p = varargin{1};
                n = varargin{2};
                a = n;
                b = bsxfun(@rdivide, bsxfun(@times, n, 1-p), p);
        end
        el = Elog(a, b);
        return
    end
    
    % Usual expected value
    % E[ln x] = psi(a) - psi(a+b)
    a  = varargin{1};
    b  = varargin{2};
    el = bsxfun(@minus, psi(a), psi(bsxfun(@plus, a, b)));
end