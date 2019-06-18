function out = Vlog(varargin)
%__________________________________________________________________________
% anatomy.math.probability.gamma.Vlog
%--------------------------------------------------------------------------
% FORMAT vl = Vlog(a, b)
% FORMAT vl = Vlog(lam, n, (K), ('normal'))
% FORMAT vl = Vlog(b,   n, a,    'gamma')
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    if ischar(varargin{end}) && ...
       startsWith(varargin{end}, 'g', 'IgnoreCase', true)
        % ----------
        % GAMMA CONJ
        % ----------
        a = bsxfun(@times, varargin{3}, varargin{2});
        
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
        
    else
        % -----
        % GAMMA
        % -----
        a = varargin{1};
    end
    
    out = psi(1,a);
end