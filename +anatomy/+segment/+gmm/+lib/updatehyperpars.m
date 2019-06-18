function [GaussPrior,extras] = updatehyperpars(cluster,GaussPrior,varargin)
%__________________________________________________________________________
% anatomy.segment.gmm.lib.updatehyperpars
%--------------------------------------------------------------------------
% FORMAT [GaussPrior,extras] = gmm.lib.updatehyperpars(cluster,GaussPrior,varargin)
%
% REQUIRED
% cluster    - 1xS cell array where cluster{s} = {{MU,b},{V,n}}
% GaussPrior - {MU0,b0,V0,n0}
%
% OPTIONAL
% constrained - Optimise hierarchical prior on V [false]
% figname     - Postfix added to figure name ['']
% verbose     - Verbosity level: [false]=quiet, true=display
%
% OUTPUT
% GaussPrior - New {MU0,b0,V0,n0}
% extras     - Struct with lower bound information, etc.
%
% Update of VB-GMM hyper-parameters (m,b,W,n).
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

import anatomy.segment.gmm.*            % plot.gaussprior
import anatomy.math.matrix.*           % posdef.logdet
import anatomy.math.probability.*      % loggamma, digamma, wishart

% Parse optional arguments
%--------------------------------------------------------------------------
p = inputParser;
p.FunctionName = 'gmm.lib.updatehyperpars';
p.addParameter('constrained',0,@islogical);
p.addParameter('figname','',@ischar);
p.addParameter('verbose',0,@islogical);
p.parse(varargin{:});
constrained = p.Results.constrained;
figname     = p.Results.figname;
verbose     = p.Results.verbose;

% Parameters
S = numel(cluster); % Number of posteriors

m0 = GaussPrior{1};
b0 = GaussPrior{2};
W0 = GaussPrior{3};
n0 = GaussPrior{4};

N = size(m0,1);
K = size(m0,2);

% pre-allocate
LogDetW0  = zeros(size(n0));
V         = zeros(size(W0));
p         = zeros(size(n0));
p0        = 0;

% -------------------------------------------------------------------------
%   Gauss-Wishart "mean" parameters
% -------------------------------------------------------------------------

for k=1:K
    
    % ---------------------------------------------------------------------
    % Update m0 (mode, closed-form)
    
    Lambda   = 0;
    LambdaMu = 0;
    for s=1:S
        [m,~,W,n] = get_posteriors(cluster,s);
        Lambda    = Lambda   + n(k)*W(:,:,k);
        LambdaMu  = LambdaMu + n(k)*W(:,:,k)*m(:,k);
    end
    m0(:,k) = Lambda \ LambdaMu;
    
    % ---------------------------------------------------------------------

    
    % ---------------------------------------------------------------------
    % Update b0 (mode, closed-form)

    b0(k)= 0;
    for s=1:S
        [m,b,W,n] = get_posteriors(cluster,s);
        m1 = m(:,k) - m0(:,k);
        b0(k) = b0(k) + m1.' * (n(k)*W(:,:,k)) * m1 + N/b(k);
    end
    b0(k) = N*S/b0(k);
    
    % ---------------------------------------------------------------------

end


% =========================================================================
% NOT CONSTRAINED
if ~constrained
    
    % ---------------------------------------------------------------------
    %   Gauss-Wishart "precision" parameters
    % ---------------------------------------------------------------------

    for k=1:K
        
        % ---
        % Set up some constants
        sumLogDet = 0;
        sumPsi    = 0;
        Wn        = 0;
        for s=1:S
            [~,~,W,n] = get_posteriors(cluster,s);
            sumLogDet = sumLogDet + posdef.logdet(W(:,:,k));
            sumPsi    = sumPsi    + digamma(n(k)/2, N);
            Wn        = Wn        + n(k)*W(:,:,k);
        end
        sumLogDet = sumLogDet/S;
        sumPsi    = sumPsi/S;
        Wn        = Wn/S;
            
        % -----------------------------------------------------------------
        % Update n0 (mode, Gauss-Newton [convex])
        E = inf;
        for gniter=1:1000
            
            % -------------------------------------------------------------
            % Update W0 (mode, closed-form)
            W0(:,:,k)   = Wn/n0(k);
            LogDetW0(k) = posdef.logdet(W0(:,:,k));
            % -------------------------------------------------------------
            
            % ---
            % Objective function
            Eprev = E;
            E = 0.5*S*n0(k)*( LogDetW0(k) - sumLogDet - sumPsi ) ...
                + S*loggamma(n0(k)/2, N);
            
            if E == Eprev
                break;
            end

            % ---
            % Gradient & Hessian
            g = 0.5*S*( LogDetW0(k) - sumLogDet - sumPsi ...
                         + digamma(n0(k)/2, N) );
            H = S/4*digamma(n0(k)/2, N, 1);

            % ---
            % Update
            n0(k) = max(n0(k) - H\g, N-1+2*eps);
            
        end
        % -----------------------------------------------------------------
    
    end
    
    % ---------------------------------------------------------------------
    %   Save results
    % ---------------------------------------------------------------------
    extras.b   = b0;
    extras.m   = m0;
    extras.n   = n0;
    extras.W   = W0;
    extras.ldW = LogDetW0;
    extras.lb  = 0;
    
    GaussPrior{1} = m0;
    GaussPrior{2} = b0;
    GaussPrior{3} = W0;
    GaussPrior{4} = n0;
    
% =========================================================================
% CONSTRAINED
else

    lb = -inf;
    for em=1:50
        % ---
        % Starting estimate
        if p0 == 0
            p0 = 0;
            V0 = 0;
            for k=1:K
                p0 = p0 + S*n0(k);
                for s=1:S
                    [~,~,W,n] = get_posteriors(cluster,s);
                    V0 = V0 + posdef.inv(n(k)*W(:,:,k));
                end
            end
            p0 = p0/K;
            V0 = V0/K;
        end

        % -----------------------------------------------------------------
        %   Gauss-Wishart "precision" parameters
        % -----------------------------------------------------------------

        for k=1:K

            % ---
            % Set up some constants
            % > compute sum E[logdet W] and sum psi(nu/2)
            logDetW  = 0;
            psiN     = 0;
            Lambda   = 0;
            for s=1:S
                [~,~,W,n] = get_posteriors(cluster,s);
                logDetW = logDetW  + posdef.logdet(W(:,:,k));
                psiN    = psiN     + digamma(n(k)/2, N);
                Lambda  = Lambda   + n(k)*W(:,:,k);
            end
            logDetW  = logDetW/S;
            psiN = psiN/S;


            % -------------------------------------------------------------
            % Update n0 (mode, Gauss-Newton [convex])
            E = inf;
            for gniter=1:100

                % ---------------------------------------------------------
                % Update {p,V} for W0 (posterior, closed form)
                p(k)       = p0 + S*n0(k);
                V(:,:,k)   = posdef.inv(posdef.inv(V0) + Lambda);
                % Useful values
                W0(:,:,k)   = posdef.inv(wishart.E(V(:,:,k), p(k)));
                LogDetW0(k) = -wishart.Elogdet(V(:,:,k), p(k));
                % ---------------------------------------------------------

                % ---
                % Objective function                
                E1 = S*n0(k)/2 * (LogDetW0(k) - logDetW - psiN) ...
                     + S*loggamma(n0(k)/2, N);
                E = [E E1];
                
                subgain = get_gain(E);                
                if subgain < 1e-6
                    % Finished
                    break
                end

                % ---
                % Gradient & Hessian
                g = S/2*(LogDetW0(k) - logDetW - psiN + digamma(n0(k)/2, N));
                H = S/4 * digamma(n0(k)/2, N, 1);

                % ---
                % Update
                n0(k) = max(n0(k) - H\g, N-1+2*eps);
            end
            % ------------------------------------------------------------

        end


        % -----------------------------------------------------------------
        %   Inverse-Wishart parameters
        % -----------------------------------------------------------------

        % ---
        % Set up some constants
        % > compute sum Logdet(psi) and sum psi(m/2)
        sumlogV = 0;
        sumPsi  = 0;
        pV      = 0;
        for k=1:K
            sumlogV = sumlogV + posdef.logdet(V(:,:,k));
            sumPsi  = sumPsi  + digamma(p(k)/2, N);
            pV      = pV      + p(k)*V(:,:,k);
        end
        sumlogV = sumlogV/K;
        sumPsi  = sumPsi/K;
        pV      = pV/K;


        % -----------------------------------------------------------------
        % Update p0 (mode, Gauss-Newton [convex])
        E = inf;
        for gniter=1:1000

            % -------------------------------------------------------------
            % Update V0 (closed-form)
            V0 = pV/p0;
            LogDetV0 = posdef.logdet(V0);
            % -------------------------------------------------------------

            % ---
            % Objective function
            Eprev = E;
            E = p0*K/2*( N*LogDetV0 - sumlogV - sumPsi ) + K*loggamma(p0/2, N);
            if E == Eprev
                break;
            end

            % ---
            % Gradient & Hessian
            g = K/2*( LogDetV0 - sumlogV - sumPsi + digamma(p0/2, N) );
            H = K/4*digamma(p0/2, N, 1);

            % ---
            % Update
            p0 = max(p0 - H\g, N-1+2*eps);

        end
        % -----------------------------------------------------------------
    
        
        % ---
        % Objective function        
        nlb  = 0;
        for k=1:K
            nlb  = nlb - wishart.kl(V(:,:,k), p(k), V0, p0);
        end
        
        lb   = [lb nlb];      
        gain = get_gain(lb);        
        
        if gain < 1e-3
            % Finished
            break
        end
        
    end % < "EM" loop

    % ---------------------------------------------------------------------
    %   Save results
    % ---------------------------------------------------------------------
    extras.b   = b0;
    extras.m   = m0;
    extras.n   = n0;
    extras.W   = W0;
    extras.ldW = LogDetW0;
    extras.V   = V;
    extras.p   = p;
    extras.V0  = V0;
    extras.p0  = p0;
    extras.lb  = 0;
    for k=1:K
        extras.lb  = extras.lb - wishart.kl(V(:,:,k), p(k), V0, p0);
    end  
    
    GaussPrior{1} = m0;
    GaussPrior{2} = b0;
    GaussPrior{3} = W0;
    GaussPrior{4} = n0;
end

if verbose
    % Visualise results
    plot.gaussprior(GaussPrior,figname);
end

% =========================================================================
function [m,b,W,n] = get_posteriors(cluster,s)
m = cluster{s}{1}{1};
b = cluster{s}{1}{2};
W = cluster{s}{2}{1};
n = cluster{s}{2}{2};

%==========================================================================
function gain = get_gain(vals)
% FORMAT gain = get_gain(vals)
vals = vals(:);
gain = abs((vals(end - 1) - vals(end))/(max(vals(isfinite(vals))) - min(vals(isfinite(vals)))));   