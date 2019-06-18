function l = kl(varargin)
%__________________________________________________________________________
% anatomy.math.probability.wishart.beta.kl
%--------------------------------------------------------------------------
% FORMAT l = kl(a1, b1, a0, b0)
% FORMAT l = kl(p1, n1, p0, n0, 'ber')
% FORMAT l = kl(p1, n1, p0, n0, k, 'bin')
% FORMAT l = kl(p1, n1, p0, n0, r, 'nbin')
% FORMAT l = kl(p1, n1, 'geom')
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    import anatomy.math.probability.beta   % kl
    
    % Check if we are in the reparameterised case
    if nargin > 4
        switch lower(varargin{end})
            case 'ber'
                p1 = varargin{1};
                n1 = varargin{2};
                p0 = varargin{3};
                n0 = varargin{4};
                a1 = bsxfun(@times, n1, p1);
                b1 = bsxfun(@times, n1, 1-p1);
                a0 = bsxfun(@times, n0, p0);
                b0 = bsxfun(@times, n0, 1-p0);
            case 'bin'
                p1 = varargin{1};
                n1 = varargin{2};
                p0 = varargin{3};
                n0 = varargin{4};
                k  = varargin{5};
                a1 = bsxfun(@times, k, bsxfun(@times, n1, p1));
                b1 = bsxfun(@times, k, bsxfun(@times, n1, 1-p1));
                a0 = bsxfun(@times, k, bsxfun(@times, n0, p0));
                b0 = bsxfun(@times, k, bsxfun(@times, n0, 1-p0));
            case 'nbin'
                p1 = varargin{1};
                n1 = varargin{2};
                p0 = varargin{3};
                n0 = varargin{4};
                r  = varargin{5};
                b1 = bsxfun(@times, r, n1);
                a1 = bsxfun(@rdivide, bsxfun(@times, p1, b1), 1-p1);
                b0 = bsxfun(@times, r, n0);
                a0 = bsxfun(@rdivide, bsxfun(@times, p0, b0), 1-p0);
            case 'geom'
                p1 = varargin{1};
                n1 = varargin{2};
                p0 = varargin{3};
                n0 = varargin{4};
                a1 = n1;
                b1 = bsxfun(@rdivide, bsxfun(@times, n1, 1-p1), p1);
                a0 = n0;
                b0 = bsxfun(@rdivide, bsxfun(@times, n0, 1-p0), p0);
        end
        l = kl(a1, b1, a0, b0);
        return
    end
    
    % Usual KL-divergence
    % KL(a1,b1||a0,b0) = ln B(a0,b0) - ln B(a1,b1) 
    %                    + (a1-a0)psi(a1) + (b1-b0)psi(b1)
    %                    + (a0-a1+b0-b1)psi(a1+b1)
    a1 = varargin{1};
    b1 = varargin{2};
    a0 = varargin{3};
    b0 = varargin{4};
    l = bsxfun(@times, ...
        bsxfun(@plus, bsxfun(@minus, a0, a1), bsxfun(@minus, b0, b1)), ...
        psi(bsxfun(@plus, a1, b1)));
    l = bsxfun(@plus,  l, bsxfun(@times, bsxfun(@minus, a1, a0), psi(a1)));
    l = bsxfun(@plus,  l, bsxfun(@times, bsxfun(@minus, b1, b0), psi(b1)));
    l = bsxfun(@plus,  l, betaln(a0,b0));
    l = bsxfun(@minus, l, betaln(a1,b1));
end