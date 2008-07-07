% PLOT - show domain on current figure
%
%  h = PLOT(d) draws many geometry features of a domain onto the current figure.
%   h is a column vector of handles to all objects plotted.
%
%  h = PLOT(d, opts) modifies this by options given by opts struct, including,
%   opts.approxp: if true, show approx polygons for each piece (default true)
%   opts.gridinside: only if >0, show gridpoints inside domain (default 0)
%
%   Also all options in SHOWSEGMENTS have effect.

function h = plot(d, o)

if nargin<2, o = []; end
if ~isfield(o, 'gridinside'), o.gridinside=0; end  % default no grid
if ~isfield(o, 'approxp'), o.approxp = 1; end  % default show polygon

h = domain.showsegments(d.seg, d.pm, o);       % show all segments

if o.approxp
  for piece=0:max(d.spiece)        % show approx poly for each connected piece
    js = find(d.spiece==piece);
    v = domain.approxpolygon(d.seg(js), d.pm(js));
    if ~isempty(v), h = [h; plot(real([v; v(1)]), imag([v; v(1)]), '--r')]; end
  end
end

l = 0.1;                           % show corner fans: radius
for j=find(~isnan(d.cloc))         % show list of valid corners only
  x = real(d.cloc(j)); y = imag(d.cloc(j));
  h = [h; plot(x, y, '.g', 'markersize', 20)];
  angfrac = 0:.05:1;
  %t = d.cangoff(j) * d.cang(j).^angfrac;      % when cang was on unit circle
  t = d.cangoff(j) * exp(1i*d.cang(j).*angfrac);
  h = [h; patch([x+l*real(t) x], [y+l*imag(t) y], 'k')]; % filled polygon patch
end

if o.gridinside>0                  % want grid of pts inside the domain
  if d.exterior
    x1 = [-2;2]; x2 = [-2;2];      % rect box to show ext domain grid
  else
    x1 =real(d.x); x2 =imag(d.x);  % x1,x2 coords of all quadr pts
  end
  [xx,yy] = meshgrid(min(x1):o.gridinside:max(x1), ...
                     min(x2):o.gridinside:max(x2));
  i = d.inside(xx+1i*yy);
  hold on; h = [h; plot(xx(i), yy(i), '.', 'markersize', 1)];
end

% show stats... the 1,1 are in case d.x=[] which happens for the entire plane
h = [h; text(max([real(d.x); 1]), max([imag(d.x); 1]), ...
             sprintf('area = %g\nperim = %g', d.area, d.perim))];
