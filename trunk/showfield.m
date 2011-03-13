% SHOWFIELD - show 2D gridded data as Re and Im images
%
%  SHOWFIELD(g, h, u, c, name) plots figure with
%   2D grid data u(:,:) with x-grid g and y-grid h. name gives the figure
%   title name, and c the size of the symmetric caxis scale. If u is complex,
%   shows Re and Im part as separate figures. If c = [], sensibel colorscales
%   are chosen
%
%  barnett 8/6/08

function showfield(g, h, u, c, name)
cfac = 1.5;
choosec = isempty(c);
re = isreal(u);                         % true if only Re needed
if re, figure('name', name); else, figure('name', sprintf('%s, Re', name)); end
imagesc(g, h, real(u)); set(gca, 'ydir', 'normal'); axis equal tight;
if choosec                          % use L4 norm to choose sensible caxis
  ug = u(find(isfinite(u)));
  c = cfac * (sum((real(ug(:)).^2).^2)/numel(ug)).^(1/4);
end
caxis(c*[-1 1]);  colormap(jet(256)); colorbar;
if ~re
  figure('name', sprintf('%s, Im', name));
  imagesc(g, h, imag(u)+0./~isnan(u)); % hack needed since imag(NaN)=0 !
  set(gca, 'ydir', 'normal'); axis equal tight;
  if choosec
    c = cfac * (sum((imag(ug(:)).^2).^2)/numel(ug)).^(1/4);
  end
  caxis(c*[-1 1]);  colormap(jet(256)); colorbar;
end
