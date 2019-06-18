function v = V(varargin)
%__________________________________________________________________________
% anatomy.math.probability.wishart.beta.V
%--------------------------------------------------------------------------
% FORMAT v = V(a, b)
% FORMAT v = V(p, n, 'ber')
% FORMAT v = V(p, n, k, 'bin')
% FORMAT v = V(p, n, r, 'nbin')
% FORMAT v = V(p, n, 'geom')
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    import anatomy.math.probability.beta   % V
    
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
        v = V(a, b);
        return
    end
    
    % Usual variance
    % V[x] = ab/[(a+b)^2*(a+b+1)]
    a  = varargin{1};
    b  = varargin{2};
    v = bsxfun(@plus, a, b);
    v = bsxfun(@times, v.^2, v+1);
    v = bsxfun(@rdivide, bsxfun(@times, a, b), v);
end