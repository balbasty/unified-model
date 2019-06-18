function l = logpdf(x, varargin)
%__________________________________________________________________________
% anatomy.math.probability.wishart.beta.logpdf
%--------------------------------------------------------------------------
% FORMAT l = logpdf(x, a, b)
% FORMAT l = logpdf(x, p, n, 'ber')
% FORMAT l = logpdf(x, p, n, k, 'bin')
% FORMAT l = logpdf(x, p, n, r, 'nbin')
% FORMAT l = logpdf(x, p, n, 'geom')
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    import anatomy.math.probability.beta   % logpdf
    
    % Check if we are in the reparameterised case
    if nargin > 3
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
        l = logpdf(x, a, b);
        return
    end
    
    % Usual log-pdf
    a   = varargin{1};
    b   = varargin{2};
    l = bsxfun(@plus,  bsxfun(@times, (a-1), log(x)), ...
                       bsxfun(@times, (b-1), log(1-x)));
    l = bsxfun(@minus, l, betaln(a, b));
end