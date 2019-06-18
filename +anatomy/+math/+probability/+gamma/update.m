function [par1, n1] = update(varargin)
%__________________________________________________________________________
% anatomy.math.probability.gamma.update
%--------------------------------------------------------------------------
% FORMAT [lam1, n1]  = update(lam, n,        lam0, n0, ('normal'))
% FORMAT [lam1, n1]  = update(ss0, ss1, ss2, lam0, n0, (mu=0), (Lambda=eye), 'normal')
% FORMAT [beta1, n1] = update(beta, n, beta0, n0, 'gamma')
% FORMAT [beta1, n1] = update(ss0, ss1,  beta0, n0, alpha, 'gamma')
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    if ischar(varargin{end}) && ...
       startsWith(varargin{end}, 'g', 'IgnoreCase', true)
        % -----
        % GAMMA
        % -----
        varargin = varargin(1:end-1);
        
        if nargin > 4
            % Sufficient statistics case
            ss0   = varargin{1};
            ss1   = varargin{2};
            beta0 = varargin{3};
            n0    = varargin{4};
            alpha = varargin{5};
            n1 = n0 + ss0;
            par1 = n0/beta0 + ss1/alpha;
            par1 = n1 / par1;
        else
            % Average case
            beta    = varargin{1};
            n       = varargin{2};
            beta0   = varargin{3};
            n0      = varargin{4};
            n1      = n + n0;
            par1 = n1 / (n0/beta0 + n/beta);
        end
        
    else
        % ------
        % NORMAL
        % ------
        if ischar(varargin{end})
            varargin = varargin(1:end-1);
        end
    
        if nargin > 4
            % Sufficient statistics case
            ss0     = varargin{1};
            ss1     = varargin{2};
            ss2     = varargin{3};
            lambda0 = varargin{4};
            n0      = varargin{5};
            K  = size(ss2, 1);
            if nargin < 7
                Lambda = eye(K);
                if nargin < 6
                    mu = zeros(size(ss1));
                else
                    mu = varargin{6};
                end
            else
                Lambda = varargin{7};
            end
            n      = ss0;
            lambda = trace(ss2*Lambda) ...
                     - 2 * mu'*Lambda*ss1 ...
                     + ss0 * mu'*Lambda*mu;
            lambda = K*ss0/lambda;
            [par1, n1] = gamma_up(lambda, n, lambda0, n0);
        else
            % Average case
            lambda  = varargin{1};
            n       = varargin{2};
            lambda0 = varargin{3};
            n0      = varargin{4};
            n1      = n + n0;
            par1 = n1 / (n0/lambda0 + n/lambda);
        end
    
    end
    
end