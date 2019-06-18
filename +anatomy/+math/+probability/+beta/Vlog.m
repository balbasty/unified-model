function vl = Vlog(varargin)
%__________________________________________________________________________
% anatomy.math.probability.wishart.beta.Vlog
%--------------------------------------------------------------------------
% FORMAT vl = Vlog(a, b)
% FORMAT vl = Vlog(p, n, 'ber')
% FORMAT vl = Vlog(p, n, k, 'bin')
% FORMAT vl = Vlog(p, n, r, 'nbin')
% FORMAT vl = Vlog(p, n, 'geom')
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    import anatomy.math.probability.beta   % Vlog
    
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
        vl = Vlog(a, b);
        return
    end
    
    % Usual variance
    % V[ln x] = psi_1(a) - psi_1(a+b)
    a  = varargin{1};
    b  = varargin{2};
    vl = bsxfun(@minus, psi(1,a), psi(1,bsxfun(@plus, a, b)));
end