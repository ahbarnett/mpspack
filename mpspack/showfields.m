% SHOWFIELDS - show multiple 2D gridded data as set of subplot images
%
%  SHOWFIELDS(g, h, u, c, name) plots figure with subplot image j from the
%   2D grid data u(:,:,j) with x-grid g and y-grid h. name gives the figure
%   title name, and c the size of the symmetric caxis scale. If u is complex,
%   shows Re and Im part separately.

function showfields(g, h, u, c, name)
n = size(u, 3);                           % # fields or subplots
re = isreal(u);                         % true if only Re needed
nh = floor(sqrt(1.8*n)); nv = ceil(n/nh); % how many across and down, subplot
if re, figure('name', name); else, figure('name', sprintf('%s, Re', name)); end
for j=1:n, subplot(nv, nh, j); imagesc(g, h, real(squeeze(u(:,:,j))));
  caxis(c*[-1 1]); set(gca, 'ydir', 'normal'); axis equal tight;
  colormap(jet(256));
end
subplotspace('vertical', -10); subplotspace('horizontal', -15);
if ~re
  figure('name', sprintf('%s, Im', name));
  for j=1:n, subplot(nv, nh, j); imagesc(g, h, imag(squeeze(u(:,:,j))));
    caxis(c*[-1 1]); set(gca, 'ydir', 'normal'); axis equal tight;
    colormap(jet(256));
  end
  subplotspace('vertical', -10); subplotspace('horizontal', -15);
end
