function p = logpdf(x, mu, varargin)
%__________________________________________________________________________
% anatomy.math.probability.normal.logpdf
%--------------------------------------------------------------------------
% FORMAT l = logpdf(x, mu,    sigma)
% FORMAT l = logpdf(x, mu,    lambda, 'precision')
% FORMAT l = logpdf(x, mu, n, sigma)
% FORMAT l = logpdf(x, mu, n, lambda, 'precision')
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    import anatomy.math.matrix
    import anatomy.math.probability.normal
    
    % Check if we are in the reparameterised case
    if nargin == 5
        if startsWith(varargin{3}, 'p', 'IgnoreCase', true)
            p = logpdf(x, mu, varargin{2}*varargin{1}, 'precision');
        else
            p = logpdf(x, mu, varargin{2}/varargin{1});
        end
        return
    elseif nargin == 4 && ~ischar(varargin{2})
        p = logpdf(x, mu, varargin{2}/varargin{1});
        return
    end
    
    % Else, set default values
    if nargin < 4
        mode = 'covariance';
        if nargin < 2
            mu = zeros(size(x));
        end
        if nargin < 3
            varargin{1} = eye(numel(mu));
        end
    else
        mode = varargin{2};
    end
    precision = startsWith(mode, 'p', 'IgnoreCase', true);
    K = size(varargin{1}, 1);
    
    if precision
        p = -0.5*( K*log(2*pi) - posdef.logdet(varargin{1}) + (x(:)-mu(:))'*varargin{1}*(x(:)-mu(:)) );
    else
        p = -0.5*( K*log(2*pi) + posdef.logdet(varargin{1}) + (x(:)-mu(:))'*(varargin{1}\(x(:)-mu(:))) );
    end
end