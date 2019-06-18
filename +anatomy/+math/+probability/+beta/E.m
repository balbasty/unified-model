function e = E(varargin)
%__________________________________________________________________________
% anatomy.math.probability.wishart.beta.E
%--------------------------------------------------------------------------
% FORMAT e = E(a, b)
% FORMAT e = E(p, n, 'ber')
% FORMAT e = E(p, n, k, 'bin')
% FORMAT e = E(p, n, r, 'nbin')
% FORMAT e = E(p, n, 'geom')
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % Check if we are in the reparameterised case
    if nargin > 2
        p = varargin{1};
        e = p;
        return
    end
    
    % Usual expected value
    % E[x] = a/(a+b)
    a = varargin{1};
    b = varargin{2};
    e = bsxfun(@rdivide, a, bsxfun(@plus, a, b));
end