function l = kl(varargin)
%__________________________________________________________________________
% anatomy.math.probability.gamma.kl
%--------------------------------------------------------------------------
% FORMAT l = kl(alpha1,  beta1, alpha0,  beta0)
% FORMAT l = kl(lambda1, n1,    lambda0, n0,    K,     ('normal'))
% FORMAT l = kl(beta1,   n1,    beta0,   n0,    alpha, 'gamma')
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % Check if we are in the reparameterised case
    if nargin == 5
        if ischar(varargin{5})
            varargin{5} = 1;
        end
        a1  = 0.5*bsxfun(@times, varargin{5}, varargin{2});
        b1  = bsxfun(@rdivide, a1, varargin{1});
        a0  = 0.5*bsxfun(@times, varargin{5}, varargin{4});
        b0  = bsxfun(@rdivide, a0, varargin{3});
    elseif nargin > 5
        if startsWith(varargin{6}, 'n', 'IgnoreCase', true)
            a1  = 0.5*bsxfun(@times, varargin{5}, varargin{2});
            b1  = bsxfun(@rdivide, a1, varargin{1});
            a0  = 0.5*bsxfun(@times, varargin{5}, varargin{4});
            b0  = bsxfun(@rdivide, a0, varargin{3});
        elseif startsWith(varargin{6}, 'g', 'IgnoreCase', true)
            a1  = bsxfun(@times, varargin{5}, varargin{2});
            b1  = bsxfun(@rdivide, a1, varargin{1});
            a0  = bsxfun(@times, varargin{5}, varargin{4});
            b0  = bsxfun(@rdivide, a0, varargin{3});
        end
    else
        a1 = varargin{1};
        b1 = varargin{2};
        a0 = varargin{3};
        b0 = varargin{4};  
    end
    
    % Usual KL
    % KL = - a0*log(b0/b1) + a1*(b0/b1 - 1) + (a1-a0)*psi(a1)
    %      + ln G(a0) - ln G(a1)
    l = bsxfun(@times, bsxfun(@minus, a1, a0), psi(a1));
    l = bsxfun(@minus, l, bsxfun(@times, a0, log(bsxfun(@rdivide, b0, b1))));
    l = bsxfun(@minus, l, bsxfun(@times, a1, bsxfun(@rdivide, b0, b1)-1));
    l = bsxfun(@plus,  l, gammaln(a0)) ;
    l = bsxfun(@minus, l, gammaln(a1)) ;
    
end