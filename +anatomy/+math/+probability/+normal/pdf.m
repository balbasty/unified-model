function p = pdf(x, mu, varargin)
%__________________________________________________________________________
% anatomy.math.probability.normal.pdf
%--------------------------------------------------------------------------
% FORMAT p = pdf(x, mu,    sigma)
% FORMAT p = pdf(x, mu,    lambda, 'precision')
% FORMAT p = pdf(x, mu, n, sigma)
% FORMAT p = pdf(x, mu, n, lambda, 'precision')
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    import anatomy.math.matrix
    import anatomy.math.probability.normal

    % Check if we are in the reparameterised case
    if nargin == 5
        if startsWith(varargin{3}, 'p', 'IgnoreCase', true)
            p = pdf(x, mu, varargin{2}.*varargin{1}, 'precision');
        else
            p = pdf(x, mu, varargin{2}./varargin{1});
        end
        return
    elseif nargin == 4 && ~ischar(varargin{2})
        p = pdf(x, mu, varargin{2}./varargin{1});
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
    
    % Usual PDF
    if precision
        p = (exp(posdef.logdet(varargin{1}/(2*pi))))^(0.5) * exp(-0.5*(x(:)-mu(:))'*varargin{1}*(x(:)-mu(:)));
    else
        p = (exp(posdef.logdet(varargin{1}*2*pi)))^(-0.5) * exp(-0.5*(x(:)-mu(:))'*(varargin{1}\(x(:)-mu(:))));
    end
    
end