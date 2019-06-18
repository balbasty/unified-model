function out = E(varargin)
%__________________________________________________________________________
% anatomy.math.probability.gamma.E
%--------------------------------------------------------------------------
% FORMAT e = E(a, b)
% FORMAT e = E(lam, n, (K), ('normal'))
% FORMAT e = E(b,   n, a,    'gamma')
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    if ischar(varargin{end}) && ...
       startsWith(varargin{end}, 'g', 'IgnoreCase', true)
        % ----------
        % GAMMA CONJ
        % ----------
            out = varargin{2};
            
    elseif ( ischar(varargin{end}) && ...
             startsWith(varargin{end}, 'n', 'IgnoreCase', true) ) || ...
           (nargin == 3)
        % -----------
        % NORMAL CONJ
        % -----------
            out = varargin{2};
            
    else
        % -----
        % GAMMA
        % -----
            out = bsxfun(@rdivide, varargin{1}, varargin{2});
    end
end