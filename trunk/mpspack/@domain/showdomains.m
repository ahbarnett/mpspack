% SHOWDOMAINS - plot all domains on current figure using a color for each
%
function h = showdomains(dlist, opts)
if nargin<2, opts = []; end
h = [];
i = 0;
for d=dlist
  hd = d.plot(opts);
  % use binary RGB sequence...
  utils.monochrome(hd, [mod(floor(i/4),2), mod(floor(i/2),2), mod(i,2)]);
  h = [h; hd];
  i = i+1;
end
