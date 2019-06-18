function p = pdf(X, varargin)
%__________________________________________________________________________
% anatomy.math.probability.wishart.pdf
%--------------------------------------------------------------------------
% FORMAT pdf = wishart_pdf(X, V,      n)
% FORMAT pdf = wishart_pdf(X, Lambda, n, 'normal')
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    import anatomy.math.probability            % loggamma
    import anatomy.math.probability.wishart    % pdf

    % Check if we are in the reparameterised case
    if nargin == 4
        p = pdf(X, varargin{1}/varargin{2}, varargin{2});
        return
    end
    
    % Usual pdf
    V   = varargin{1};
    n   = varargin{2};
    K   = size(V, 1);
    p = det(X)^(n-K-1) * exp(-0.5*trace(V\X)) ...
          / ( 2^(n*K/2) * det(V)^(n/2) * exp(loggamma(0.5*n, K)) );
    
end