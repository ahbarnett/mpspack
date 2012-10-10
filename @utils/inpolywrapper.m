function i = inpolywrapper(p, v)
% INPOLYWRAPPER - logical array of if numbers p are in polygon v (various codes)
%
% i = INPOLYWRAPPER(p, v) where p and v are column vectors of complex numbers,
%   returns logical array i of size same as p, giving whether each point p is
%   within simply-connected polygon with vertices v.
%
% Chooses between one of 3 methods:
%   * Matlab's native inpolygon
%        inpolygon_vec is very slow, 1.5e7 point-vertices per sec on i7-3720QM
%        But it does seem to do an initial bounding-box test.
%   * Darren Engwirda's inpoly
%        Claims to have some intelligent sorting of points in p.
%        Around 10-20x faster than matlab native.
%   * MEX interface to Wm. Randolph Franklin's C code
%        Brute forces it and loops over points in p. Around 30x matlab native.
%   * MEX interface to Bruno Luong's C code
%        Has bunch of sorting stuff. Excellent: 1.5e9 pt-verts/sec on i7-3720QM
%        (and eg 2e10 if only 10% lie in bb!). Ie 100x matlab native.
%
% See also: INPOLYGON, TEST/TESTINPOLYWRAPPER, UTILS/INPOLYC

% Copyright (C) 2008 - 2012, Timo Betcke, Alex Barnett.
% With code by Peter Simon, 10/4/2012 to use Bruno Luong's
%   insidepoly MEX implementation

% trivial cases
if isempty(v) | numel(v)==1, i = 0*p; return; end   % otherwise bb=[], crashes
if isempty(p), i = []; return; end

% use crude initial bounding-box test for speed.... helps all the C codes!
bb = [min(real(v)) max(real(v)) min(imag(v)) max(imag(v))]; % v's bounding box
i = (real(p)>=bb(1) & real(p)<=bb(2) & imag(p)>=bb(3) & imag(p)<=bb(4));
%i = true(size(p)); % ...or bypass the bb test

% Please uncomment one of the following 4 calls... (note p(i) is pts in bb)

%   Matlab's native inpolygon (use if no MEX files work)
%i(i) = inpolygon(real(p(i)), imag(p(i)), real(v), imag(v));

%   Darren Engwirda's inpoly
%i(i) = utils.inpoly([real(p(i(:))),imag(p(i(:)))], [real(v),imag(v)]);

%   MEX interface to Wm. Randolph Franklin's C code:
%i(i) = logical(utils.inpolyc(p(i), v));   % converts int to logical

%   MEX interface to Bruno Luong's C code: (recommended)
i(i) = utils.insidepoly(real(p(i)), imag(p(i)), real(v), imag(v));
