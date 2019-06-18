function out = V(varargin)
%__________________________________________________________________________
% anatomy.math.probability.gamma.V
%--------------------------------------------------------------------------
% FORMAT v = V(a, b)
% FORMAT v = V(lam, n, (K), ('normal'))
% FORMAT v = V(b,   n, a,    'gamma')
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    if ischar(varargin{end}) && ...
       startsWith(varargin{end}, 'g', 'IgnoreCase', true)
        % ----------
        % GAMMA CONJ
        % ----------
        a = bsxfun(@times, varargin{3}, varargin{2});
        b = bsxfun(@rdivide, a, varargin{1});
                             
    elseif ( ischar(varargin{end}) && ...
             startsWith(varargin{end}, 'n', 'IgnoreCase', true) ) || ...
           (nargin == 3)
        % -----------
        % NORMAL CONJ
        % -----------
        if ischar(varargin{3})
            varargin{3} = 1;
        end
        a = 0.5 * bsxfun(@times, varargin{3}, varargin{2});
        b = bsxfun(@rdivide, a, varargin{1});
                                 
    else
        % -----
        % GAMMA
        % -----
        a = varargin{1};
        b = varargin{2};
        
    end
    out = bsxfun(@rdivide, a, b.^2);
end