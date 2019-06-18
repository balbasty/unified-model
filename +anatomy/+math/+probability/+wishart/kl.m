function l = kl(varargin)
%__________________________________________________________________________
% anatomy.math.probability.wishart.kl
%--------------------------------------------------------------------------
% FORMAT kl = kl(V1,      n1, V0,      n0)
% FORMAT kl = kl(lambda1, n1, lambda0, n0, 'normal')
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    import anatomy.math.matrix                 % posdef.logdet
    import anatomy.math.probability            % loggamma, digamma
    import anatomy.math.probability.wishart    % logpdf
    
    % Check if we are in the reparameterised case
    if nargin == 5
        l = kl(varargin{1}/varargin{2}, varargin{2}, ...
               varargin{3}/varargin{4}, varargin{4});
        return
    end
    
    % Usual KL
    V1 = varargin{1};
    n1 = varargin{2};
    V0 = varargin{3};
    n0 = varargin{4};
    K  = size(V1, 1);
    l =   0.5*n0*(posdef.logdet(V0) - posdef.logdet(V1)) ...
         + 0.5*n1*(trace(V0\V1) - K) ...
         + 0.5*(n1 - n0)*digamma(0.5*n1, K) ...
         + loggamma(0.5*n0, K) - loggamma(0.5*n1, K);
    
end