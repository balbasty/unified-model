function [p1,n1] = update(varargin)
%__________________________________________________________________________
% anatomy.math.probability.wishart.beta.update
%--------------------------------------------------------------------------
% FORMAT [p1,n1] = update(p, n, p0, n0)
%
% FORMAT [p1,n1] = update(s0, s1, p0, n0,    'suffstat', 'ber')   p=s1/s0
% FORMAT [p1,n1] = update(s0, s1, p0, n0, k, 'suffstat', 'bin')   p=s1/(k*s0)
% FORMAT [p1,n1] = update(s0, s1, p0, n0, r, 'suffstat', 'nbin')  p=s1/(r*s0+s1)
% FORMAT [p1,n1] = update(s0, s1, p0, n0,    'suffstat', 'geom')  p=s0/s1
% FORMAT [p1,n1] = update(s0, s1, p0, n0,    'suffstat', 'fgeom') p=s0/(s1+s0)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging


    if ischar(varargin{end})
        switch lower(varargin{end})
            case 'ber'
                s0 = varargin{1};
                s1 = varargin{2};
                n = s0;
                p = s1/s0;
            case 'bin'
                s0 = varargin{1};
                s1 = varargin{2};
                k  = varargin{5};
                n = s0;
                p = s1/(k*s0);
            case 'nbin'
                s0 = varargin{1};
                s1 = varargin{2};
                r  = varargin{5};
                n = s0;
                p = s1/(r*s0+s1);
            case 'geom'
                s0 = varargin{1};
                s1 = varargin{2};
                n = s0;
                p = s0/s1;
            case 'fgeom'
                s0 = varargin{1};
                s1 = varargin{2};
                n = s0;
                p = s0/(s1+s0);
        end
    else
        p = varargin{1};
        n = varargin{2};
    end
    p0 = varargin{3};
    n0 = varargin{4};
    n1 = n0+n;
    p1 = (n0*p0+n*p)/n1;

end