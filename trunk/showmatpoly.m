function showmatpoly(M, name)
% SHOWMATPOLY - show subplots containing matrix polynomial coeffs & alpha powers
%
%   Currently all plots are set to the same, symmetric, color axis, and only
%   real part is shown. Precisely-zero values are highlighted as NaNs (care!)
c = max(abs(M(:)));
g = gcf; figure(g); if nargin>1, set(g, 'name', name); end
P = size(M,3);              % polynomial order + 1
ind3 = ceil((P-1)/2)+1;     % formula from basis/evalunitcellcopies for offset
for i=1:P
  tsubplot(1,P,i);
  Mkillzeros = zeros(size(M(:,:,i))); Mkillzeros(find(M(:,:,i)==0))=NaN;
  imagesc(real(M(:,:,i))+Mkillzeros); caxis(c*[-1 1]); colormap(jet(256));
  if i==P, title(sprintf('\\alpha^{%d} caxis=%.2g', i-ind3, c)); else
  title(sprintf('\\alpha^{%d}', i-ind3)); end
end
