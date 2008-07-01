% SHOWSEGMENTS - helper which plots signed segment list
%
function h = showsegments(segs, pm)

for j=1:numel(segs)
  s = segs(j);
  h(j) = plot(s, pm(j));
  v = domain.approxpolygon(segs, pm);
  plot([v; v(1)], '--r');
end






