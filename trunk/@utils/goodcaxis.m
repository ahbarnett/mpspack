function c = goodcaxis(u)
% GOODCAXIS - use quantile of array to choose a good symm caxis for wave images
%
% Based on idea from code imgs_smart_caxis of B. Gustavsson 2005-02-09

% Copyright (C) 2010 - 2011, Alex Barnett
alpha = 0.95;                             % upper quantile (out of 1)
v = abs(u(find(~isnan(u)))); v = v(:);
[b x] = hist(v, unique(v));
ch = cumsum(b)/numel(v);
ic = find(ch > alpha);
if isempty(ic), caxis([-1 1]); warning('entire array has same value'); % give up
else, c = x(ic(1));
  caxis(c*[-1 1]);
end
