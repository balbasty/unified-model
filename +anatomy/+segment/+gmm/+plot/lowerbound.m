function lowerbound(lb, figname)
%__________________________________________________________________________
% anatomy.segment.gmm.plot.lowerbound
%--------------------------------------------------------------------------
% gmm.plot.lowerbound(lb, (wintitle))
% 
% Plot lower bound
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging


% -------------------------------------------------------------------------
% Get figure (create if it does not exist)
if nargin < 2
    figname = 'Plot GMM Lower Bound';
end
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   
clf(f);

% -------------------------------------------------------------------------
% Choose type
if isfield(lb, 'B')
    nrow = 2;
    ncol = 4;
else
    nrow = 2;
    ncol = 3;
end

% -------------------------------------------------------------------------
% Plots
subplot(nrow, ncol, sub2ind([ncol nrow], 1, 1));
plot(lb.sum)
title('Lower Bound')
subplot(nrow, ncol, sub2ind([ncol nrow], 2, 1));
if isfield(lb, 'B')
    plot(sum(lb.X,1) + sum(lb.XB,1));
else
    plot(sum(lb.X,1))
end
box on
title('Observations (E)')
subplot(nrow, ncol, sub2ind([ncol nrow], 3, 1));
plot(sum(lb.Z,1))
box on
title('Responsibilities (KL)')
subplot(nrow, ncol, sub2ind([ncol nrow], 1, 2));
plot(sum(lb.P,1))
box on
title('Proportions (KL)')
subplot(nrow, ncol, sub2ind([ncol nrow], 2, 2));
plot(sum(lb.MU,1))
box on
title('Means (KL)')
subplot(nrow, ncol, sub2ind([ncol nrow], 3, 2));
plot(sum(lb.A,1))
box on
title('Precisions (KL)')
if isfield(lb, 'B')
    subplot(nrow, ncol, sub2ind([ncol nrow], 4, 2));
    plot(sum(lb.B,1))
    box on
    title('Bias Prior')
end
drawnow