function i = inpolywrapper(p, v)
% INPOLYWRAPPER - logical array of if numbers p are in polygon v (various codes)
%
% i = INPOLYWRAPPER(p, v) where p and v are column vectors of complex numbers,
%   returns logical array i of size same as p, giving whether each point p is
%   within simply-connected polygon with vertices v.
%
% By default, chooses between one of 2 methods: a very fast MEX interface
%  due to Bruno Luong (around 1e9 pts-verts/sec), or if that's unavailable,
%  MATLAB's native inpolygon (1e7 pts-verts/sec, even on R2015).
%
% In fact, a total of 4 inpoly methods are kept, for completeness:
%   * Matlab's native inpolygon
%        inpolygon_vec is very slow, 1.5e7 point-vertices per sec on i7-3720QM
%        But it does seem to do an initial bounding-box test. R2015b no faster
%        than R2012a!
%   * Darren Engwirda's inpoly
%        Claims to have some intelligent sorting of points in p.
%        Around 10-20x faster than matlab native.
%   * MEX interface to Wm. Randolph Franklin's C code
%        Brute forces it and loops over points in p. Around 30x matlab native.
%   * MEX interface to Bruno Luong's M-file and shipped .mex object files.
%        Has bunch of sorting stuff. Blinding: 1.5e9 pt-verts/sec on i7-3720QM
%        (and eg 2e10 if only 10% lie in bb). Ie 100x matlab native!
%
% See also: INPOLYGON, TEST/TESTINPOLYWRAPPER, UTILS/INPOLYC

% Copyright (C) 2008 - 2016, Alex Barnett.
% With code by Peter Simon, 10/4/2012 to use Bruno Luong's
%   insidepoly MEX implementation
% Barnett changed to MATLAB's native inpolgon by default 10/24/13.
% switch using MEX-file check 4/12/16

% trivial cases
if isempty(v) | numel(v)==1, i = 0*p; return; end   % otherwise bb=[], crashes
if isempty(p), i = []; return; end

% use crude initial bounding-box test for speed.... helps all the C codes!
bb = [min(real(v)) max(real(v)) min(imag(v)) max(imag(v))]; % v's bounding box
i = (real(p)>=bb(1) & real(p)<=bb(2) & imag(p)>=bb(3) & imag(p)<=bb(4));
%i = true(size(p)); % ...or bypass the bb test

% Please uncomment one of the following 4 calls... (note p(i) is pts in bb)
% I recommend either the first (works out of the box), or last (for speed).
% To use the last two, you need to compile via make in @utils.

if exist('@utils/insidepoly_dblengine')==3        % MEX avail
  %   Bruno Luong's M-file and associated MEX engines: (recommended)
  i(i) = utils.insidepoly(real(p(i)), imag(p(i)), real(v), imag(v));
else
  %   MATLAB's native inpolygon (fall-back; use if no MEX files work)
  i(i) = inpolygon(real(p(i)), imag(p(i)), real(v), imag(v));
end

% for posterity only:

%   Darren Engwirda's inpoly (MATLAB)
%i(i) = utils.inpoly([real(p(i(:))),imag(p(i(:)))], [real(v),imag(v)]);

%   MEX interface to Wm. Randolph Franklin's C code:
%i(i) = logical(utils.inpolyc(p(i), v));   % converts int to logical
