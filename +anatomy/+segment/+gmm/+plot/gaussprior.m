function gaussprior(GaussPrior,figname)
%__________________________________________________________________________
% anatomy.segment.gmm.plot.gaussprior
%--------------------------------------------------------------------------
% gmm.plot.gaussprior(GaussPrior, (wintitle))
% 
% Plot VB-GMM hyper-parameters
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

figname0 = 'GaussPrior';
if nargin==2
    figname0 = [figname0 ' ' figname];
end

% -------------------------------------------------------------------------
% Get figure (create if it does not exist)
f = findobj('Type', 'Figure', 'Name', figname0);
if isempty(f)
    f = figure('Name', figname0, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   
clf(f);

m0 = GaussPrior{1};
W0 = GaussPrior{3};
n0 = GaussPrior{4};
    
P      = size(m0,1);
K      = size(m0,2);
colors = hsv(K);

MU = m0;
A  = bsxfun(@times,W0,reshape(n0,[1 1 K]));

% -------------------------------------------------------------------------
% For each input dimension
for p=1:P
    % ---------------------------------------------------------------------
    % Plot histogram and marginal density
    if P>1
        subplot(2, P, p)
    end
    hold on

    xlims = [inf -inf];
    % -----------
    % GMM Density
    for k=1:K
        x = linspace(MU(p,k)-3*A(p,p,k)^(-0.5),MU(p,k)+3*A(p,p,k)^(-0.5),100);
        y = 1/K*normpdf(x, MU(p,k), A(p,p,k)^(-0.5));
        plot(x, y, 'Color', colors(k,:), 'LineWidth', 1)
        xlims = [min([xlims(1) x]) max([xlims(2) x])];
    end
    xlabel(sprintf('x%d',p))
    ylabel('density')
    xlim(xlims);
    box on
    hold off

    % ---------------------------------------------------------------------
    % Plot joint density (X1 vs Xj)
    if p > 1
        subplot(2, P, P+p)
        hold on
        for k=1:K
            Mu1     = MU([1 p],k);
            Sigma2  = spm_matcomp('Inv', A([1 p],[1 p],k));
            Sigma   = sqrt(Sigma2);
            [x1,x2] = meshgrid(linspace(Mu1(1)-3*Sigma(1,1),Mu1(1)+3*Sigma(1,1),100)', ...
                               linspace(Mu1(2)-3*Sigma(2,2),Mu1(2)+3*Sigma(2,2),100)');
            y = mvnpdf([x1(:) x2(:)],Mu1',Sigma2);
            contour(x2, x1, reshape(y, [100 100])', 1, 'color', colors(k,:), 'LineWidth', 1);
        end
        xlabel(sprintf('x%d',p))
        ylabel('x1')
        xlim(xlims); % use same scale as histogram plot for comparison
        box on
        hold off
    end
end

drawnow