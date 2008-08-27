function showmatpoly(M, name)
% SHOWMATPOLY - show five subplots containing matrix polynomial from -2,-1,..,2
if size(M,3)~=5, error('must be a MxNx5 3d array!'); end
c = max(abs(M(:)));
g = gcf; figure(g); if nargin>1, set(g, 'name', name); end
for i=1:5
  tsubplot(1,5,i);
  imagesc(real(squeeze(M(:,:,i)))); caxis(c*[-1 1]); colormap(jet(256));
  if i==5, title(sprintf('\\alpha^{%d} caxis=%.2g', i-3, c)); else
  title(sprintf('\\alpha^{%d}', i-3)); end
end
