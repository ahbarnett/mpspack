% PLOT - plots a directed segment on current figure, using its quadrature pts
%
%  h = PLOT(seg) plots a segment, quadrature points, etc. Also plots an array
%   of segment handles with their natural senses, and number labels.
%
%  h = PLOT(seg, pm) plots segment if pm=1, or its reversal if pm=-1.
%
%  h = PLOT(s, pm, opts) allows various options in opts, namely
%   opts.arrow: if true, show direction via an arrow (default true)
%   opts.normals: if true, show directions of the normals (default true)
%   opts.blobs: if true, show quadrature pt blobs (default true)
%
% See also: pointset/PLOT, domain/SHOWSEGMENTS

function h = plot(s, pm, o)
if nargin<2, pm = 1; end                       % default sense is positive
if nargin<3, o = []; end
if ~isfield(o, 'arrow'), o.arrow = 1; end % default is show arrow
if ~isfield(o,'normals'), o.normals=1; end % default is show normals
lt = '.-';                                 % default is show quad pt blobs
if isfield(o,'blobs') & ~o.blobs, lt='-'; end % switch off blobs
closed = (abs(s.Z(0)-s.Z(1))<1e-15);       % hack to tell if segment is closed

g = gcf;
figure(g); hold on;

if numel(s)>1                       % vectorize using domain routine
  domain.showsegments(s, pm, o);
else                                % just one seg, plot it!
  if closed, h = plot(real([s.x; s.x(1)]), imag([s.x; s.x(1)]), lt);
  else, h = plot(real(s.x), imag(s.x), lt);
  end
  if o.normals,
    l = 0.1;                                       % show normals... length
    h = [h; plot([s.x(:).'; (s.x(:)+l*pm*s.nx(:)).'], 'k-')]; % uses sign from pm
  end
    
  if o.arrow
    t = 0.5 + pm * 0.04 * [-1+1i, 1, -1-1i];     % t values of arrow, with direc
    h = [h; plot(s.Z(t), '-')];
  end
end
axis equal;
hold off;
