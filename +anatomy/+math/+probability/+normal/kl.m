function l = kl(mu1, par1, mu0, par0, varargin)
%__________________________________________________________________________
% anatomy.math.probability.normal.kl
%--------------------------------------------------------------------------
% FORMAT l = kl(mu1, sigma1,  mu0, sigma0)
% FORMAT l = kl(mu1, lambda1, mu0, lambda0, 'precision')
% FORMAT l = kl(mu1, n1,      mu0, n0,      sigma)
% FORMAT l = kl(mu1, n1,      mu0, n0,      lambda, 'precision')
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    import anatomy.math.matrix
    import anatomy.math.probability.normal
    
    % Check if we are in the reparameterised case
    if nargin == 6
        K = size(varargin{1}, 1);
        if startsWith(varargin{2}, 'p', 'IgnoreCase', true)
            l = 0.5 * ( K * (par0/par1 - log(par0/par1) - 1) ...
                         + par0 * (mu0(:)-mu1(:))'*varargin{1}*(mu0(:)-mu1(:)) );
        else
            l = 0.5 * ( K * (par0/par1 - log(par0/par1) - 1) ...
                         + par0 * (mu0(:)-mu1(:))'*(varargin{1}\(mu0(:)-mu1(:))) );
        end
        return
    elseif nargin == 5 && ~ischar(varargin{1})
        K  = size(varargin{1}, 1);
        l = 0.5 * ( K * (par0/par1 - log(par0/par1) - 1) ...
                     + par0 * (mu0(:)-mu1(:))'*(varargin{1}\(mu0(:)-mu1(:))) );
        return
    end
    
    % Else, set default values
    if nargin < 5
        mode = 'covariance';
    else
        mode = varargin{1};
    end
    precision = startsWith(mode, 'p', 'IgnoreCase', true);
    K = size(par1, 1);
    
    % Common KL-divergence
    if precision
        l = 0.5 * ( trace(par1\par0) ...
                     - posdef.logdet(par0) ...
                     + posdef.logdet(par1) ...
                     - K ...
                     + (mu0(:)-mu1(:))'*par0*(mu0(:)-mu1(:)) );
    else
        l = 0.5 * ( trace(par0\par1) ...
                     - posdef.logdet(par1) ...
                     + posdef.logdet(par0) ...
                     - K ...
                     + (mu0(:)-mu1(:))'*(par1\(mu0(:)-mu1(:))) );
    end
    
end