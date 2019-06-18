function categorical(Z,Template,figname)
%__________________________________________________________________________
% anatomy.segment.gmm.plot.categorical
%--------------------------------------------------------------------------
% gmm.plot.categorical(dm, Z, Template, (wintitle))
%
% Plot (categorical) responsibilities and template (if available)
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

% -------------------------------------------------------------------------
% Get figure (create if it does not exist)
if nargin<4
    figname = 'Plot GMM Categorical';
end
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   
clf(f);

K        = size(Z,4);
colors   = hsv(K);
do_subpl = isequal(size(Z),size(Template));

% -------------------------------------------------------------------------
% Plots
if do_subpl
    subplot(121)
end

c = math.gmm.lib('plot', 'cat2rgb', Z, colors);
c = squeeze(c(:,:,:,:));
c = permute(c,[2 1 3]);
imagesc(c); axis image off xy;   
title('Z')

colormap(colors);
cb = colorbar;
set(gca, 'clim', [0.5 K+0.5]);
set(cb, 'ticks', 1:K, 'ticklabels', {1:K}); 

if do_subpl
    subplot(122)    
    
    c = cat2rgb(Template, colors);
    c = squeeze(c(:,:,:,:));
    c = permute(c,[2 1 3]);
    imagesc(c); axis image off xy;   
    
    title('Template')
    
    colormap(colors);
    cb = colorbar;
    set(gca, 'clim', [0.5 K+0.5]);
    set(cb, 'ticks', 1:K, 'ticklabels', {1:K});     
end

drawnow

% =========================================================================
function c = cat2rgb(f, pal)
% FORMAT c = cat2rgb(f, pal)
% f   - categorical (4D) image.
% pal - palette (Mx3 array or handle to palette function) [hsv]
%
% Generate an RGB volume from a categorical (e.g. responsibilities) volume.

if nargin < 2
    pal = @hsv;
end

if size(f,3)>1
    z = floor(size(f,3)/2) + 1;
    f = f(:,:,z,:);
end

tri = false;
if numel(size(f)) == 4 && size(f, 3) == 1
    tri = true;
    dm  = [size(f) 1 1];
    f   = reshape(f, [dm(1:2) dm(4)]);
end
if isa(pal, 'function_handle')
    pal = pal(size(f,3));
end

dm = [size(f) 1 1];
c   = zeros([dm(1:2) 3]); % output RGB image
s   = zeros(dm(1:2));     % normalising term

for k=1:dm(3)
    s = s + f(:,:,k);
    color = reshape(pal(k,:), [1 1 3]);
    c = c + bsxfun(@times, f(:,:,k), color);
end
if dm(3) == 1
    c = c / max(1, max(s(:)));
else
    c = bsxfun(@rdivide, c, s);
end

if tri
    c = reshape(c, [size(c, 1) size(c, 2) 1 size(c, 3)]);
end