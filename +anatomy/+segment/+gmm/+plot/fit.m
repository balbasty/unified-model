function fit(obs, cluster, PI, figname)
%__________________________________________________________________________
% anatomy.segment.gmm.plot.fit
%--------------------------------------------------------------------------
% gmm.plot.fit({X,W}, {MU,A}, PI, (wintitle))
% 
% Plot mixture fit
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

X  = [];
W  = [];
MU = [];
A  = [];

if ~iscell(obs)
    X = obs;
else
    if numel(obs) >= 1
        X = obs{1};
        if numel(obs) >= 2
            W = obs{2};
        end
    end
end
if ~iscell(cluster)
    MU = cluster;
else
    if numel(cluster) >= 1
        MU = cluster{1};
        if numel(cluster) >= 2
            A = cluster{2};
        end
    end
end

% -------------------------------------------------------------------------
% Get figure (create if it does not exist)
if nargin < 4
    figname = 'Plot GMM';
end
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   
clf(f);

% -------------------------------------------------------------------------
% Sizes / colors
P = size(X, 2);
K = size(MU, 2);
colors = hsv(K);

% -------------------------------------------------------------------------
% For each input dimension
for p=1:P
    % ---------------------------------------------------------------------
    % Plot histogram and marginal density
    subplot(2, P, p)
    hold on
    % ---------
    % Histogram
    if any(W > 1)
        % Input is already an histogram
        X1 = X(isfinite(X(:,p)),p);
        weights = W(isfinite(X(:,p)));
        [idx,edges] = discretize(X1, 64);
        centres = (edges(2:end) + edges(1:end-1))/2;
        idx2 = zeros(numel(idx), max(idx));
        idx2(sub2ind(size(idx2), 1:numel(idx), idx')) = 1;
        weights = sum(bsxfun(@times, weights, idx2), 1);
        clear idx1 idx2
        minx  = min(X1);
        maxx  = max(X1);
        weights = weights ./ (sum(weights)*(maxx-minx)/numel(weights));
        bar(centres, weights, 'EdgeColor', 'none', 'FaceColor', [0.7 0.7 0.7]);
        ymax = max(weights);
    else
        % Input is a list of observations
        [H, edges] = histcounts(X(:,p), 64, 'Normalization', 'pdf');
        centres = (edges(1:end-1) + edges(2:end))/2;
        bar(centres, H, 'EdgeColor', 'none', 'FaceColor', [0.7 0.7 0.7]);
        ymax = max(H);
    end
    xlims = [inf -inf];
    % -----------
    % GMM Density
    for k=1:K
        x = linspace(MU(p,k)-3*A(p,p,k)^(-0.5),MU(p,k)+3*A(p,p,k)^(-0.5),100);
        y = PI(k)*normpdf(x, MU(p,k), A(p,p,k)^(-0.5));
        plot(x, y, 'Color', colors(k,:), 'LineWidth', 1)
        xlims = [min([xlims(1) x]) max([xlims(2) x])];
    end
    xlabel(sprintf('x%d',p))
    ylabel('density')
    xlim(xlims);
    ylim([0 ymax]);
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

% -------------------------------------------------------------------------
% Plot proportions in the remaining square
subplot(2, P, P+1)
hold on
for k=1:K
    bar(k, PI(k), 'FaceColor', colors(k,:));
end
xlabel('class')
ylabel('proportion')
box on
hold off
drawnow