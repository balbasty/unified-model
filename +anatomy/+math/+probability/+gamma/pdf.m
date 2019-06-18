function p = pdf(x, varargin)
%__________________________________________________________________________
% anatomy.math.probability.gamma.pdf
%--------------------------------------------------------------------------
% FORMAT p = pdf(x, alpha,  beta)
% FORMAT p = pdf(x, lambda, n,    K,     ('normal'))
% FORMAT p = pdf(x, beta,   n,    alpha, 'gamma')
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % Check if we are in the reparameterised case
    if nargin == 4
        if ischar(varargin{3})
            varargin{3} = 1;
        end
        a   = 0.5 * bsxfun(@times, varargin{3}, varargin{2});
        b   = bsxfun(@rdivide, a, varargin{1});
    elseif nargin > 4 
        if startsWith(varargin{4}, 'n', 'IgnoreCase', true)
            a   = 0.5 * bsxfun(@times, varargin{3}, varargin{2});
            b   = bsxfun(@rdivide, a, varargin{1});
        elseif startsWith(varargin{4}, 'g', 'IgnoreCase', true)
            a   = bsxfun(@times, varargin{3}, varargin{2});
            b   = bsxfun(@rdivide, a, varargin{1});
        end
    else
        a = varargin{1};
        b = varargin{2};
    end
    
    % Usual pdf
    p = gampdf(x, a, 1./b);
    
end