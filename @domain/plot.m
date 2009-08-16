function h = plot(d, o)
% PLOT - show domain (or list of domains) on current figure
%
%  h = PLOT(d) draws many geometry features of a domain onto the current figure.
%   h is a column vector of handles to all objects plotted.
%
%  h = PLOT(d, opts) modifies this by options given by opts struct, including,
%   opts.approxp: if true, show approx polygons for each piece (default true)
%   opts.gridinside: only if >0, show gridpoints inside domain (default 0)
%
%   Also all options in SHOWSEGMENTS have effect.
%
% Also see: DOMAIN.SHOWDOMAINS which is the correct code for domain lists

% Copyright (C) 2008, 2009, Alex Barnett, Timo Betcke

if nargin<2, o = []; end
if numel(d)>1, domain.showdomains(d, o); return; end
% the rest of code handles a single domain object...

if ~isfield(o, 'gridinside'), o.gridinside=0; end  % default no grid
if ~isfield(o, 'approxp'), o.approxp = 1; end  % default show polygon
if ~isfield(o, 'filled'), o.filled=0; end   % to show domains as solid ?? to do

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

if o.gridinside>0                  % want grid of pts inside the domain?
  [zz] = d.grid(o.gridinside);
  h = [h; plot(real(zz), imag(zz), '.', 'markersize', 1)];
end

% show stats... the 1,1 are in case d.x=[] which happens for the entire plane
h = [h; text(max([real(d.x); 1]), max([imag(d.x); 1]), ...
             sprintf('area = %g\nperim = %g', d.area, d.perim))];

if d.refr_ind ~= 1.0
  h = [h; text(real(d.center), imag(d.center), sprintf('n=%.3g', d.refr_ind))];
end
