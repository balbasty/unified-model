function l = logpdf(X, varargin)
%__________________________________________________________________________
% anatomy.math.probability.wishart.logpdf
%--------------------------------------------------------------------------
% FORMAT l = logpdf(X, V,      n)
% FORMAT l = logpdf(X, Lambda, n, 'normal')
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    import anatomy.math.matrix                 % posdef.logdet
    import anatomy.math.probability            % loggamma
    import anatomy.math.probability.wishart    % logpdf
    
    % Check if we are in the reparameterised case
    if nargin == 4
        l = logpdf(X, varargin{1}/varargin{2}, varargin{2});
        return
    end
    
    % Usual pdf
    V   = varargin{1};
    n   = varargin{2};
    K   = size(V, 1);
    l =   0.5*(n-K-1)*posdef.logdet(X) ...
          - 0.5*trace(V\X) ...
          - 0.5*n*K*log(2) ...
          - 0.5*n*posdef.logdet(V) ...
          - loggamma(0.5*n, K);
    
end