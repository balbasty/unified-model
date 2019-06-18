function p = pdf(x, varargin)
%__________________________________________________________________________
% anatomy.math.probability.wishart.beta.pdf
%--------------------------------------------------------------------------
% FORMAT p = pdf(x, a, b)
% FORMAT p = pdf(x, p, n, 'ber')
% FORMAT p = pdf(x, p, n, k, 'bin')
% FORMAT p = pdf(x, p, n, r, 'nbin')
% FORMAT p = pdf(x, p, n, 'geom')
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    import anatomy.math.probability.beta   % pdf

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
        p = pdf(x, a, b);
        return
    end
    
    % Usual pdf
    a   = varargin{1};
    b   = varargin{2};
    p = betapdf(x, a, b);
end