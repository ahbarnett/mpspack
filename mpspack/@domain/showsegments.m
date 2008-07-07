% SHOWSEGMENTS - plot signed segment list to current figure (domain helper)
%
%  h = SHOWSEGMENTS(s, pm) plots segment pointers s with senses pm to
%   current figure. h is a column vector of handles to all objects plotted.
%
%  h = SHOWSEGMENTS(s, pm, opts) allows options in opts, namely
%   opts.arrow: if true, put an arrow on each segment (default true)
%   opts.label: if true, label each seg with its number in list (default true)
%
%  See also: domain/PLOT, segment/PLOT

function h = showsegments(segs, pm, o)
if nargin<3, o = []; end
if ~isfield(o, 'label'), o.label = 1; end      % default is label

h = [];
for j=1:numel(segs)
  s = segs(j);
  h = [h; plot(s, pm(j), o)]; hold on;          % pass opts to segment.plot
  if o.label
    xm = s.Z(1/2);                   % label at location half way along seg
    h = [h; text(real(xm), imag(xm), sprintf('%d',j))];
  end
end
