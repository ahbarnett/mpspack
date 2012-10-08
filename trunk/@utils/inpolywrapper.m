function i = inpolywrapper(p, v)
% INPOLYWRAPPER - logical array of if numbers p are in polygon v (various codes)
%
% i = INPOLYWRAPPER(p, v) where p and v are column vectors of complex numbers,
%   returns logical array i of size same as p, giving whether each point p is
%   within simply-connected polygon with vertices v.
%
% Chooses between one of 3 methods:
%   * Matlab's native inpolygon
%        inpolygon_vec is very slow, 0.5us per point per vertex (2GHz Core Duo)
%   * Darren Engwirda's inpoly
%        claims to have some intelligent sorting of points in p, but I find
%        only 1.25 times speed of matlab. Or, now, 10x slower! (2012, i7)
%   * MEX interface to Wm. Randolph Franklin's C code
%        brute forces it and loops over points in p, seems to be 70 times
%        faster than matlab, ie only 8 ns per point per vertex! (2GHz Core Duo)
%
% See also: INPOLYGON, TEST/TESTINPOLYWRAPPER, UTILS/INPOLYC

% Copyright (C) 2008 - 2012, Timo Betcke, Alex Barnett

% use crude initial bounding-box test for speed....
bb = [min(real(v)) max(real(v)) min(imag(v)) max(imag(v))]; % v's bounding box
i = (real(p)>=bb(1) & real(p)<=bb(2) & imag(p)>=bb(3) & imag(p)<=bb(4));

% uncomment one of the following calls... (note p(i) is pts in bb)

%   Matlab's native inpolygon
%i(i) = inpolygon(real(p(i)), imag(p(i)), real(v), imag(v));

%   Darren Engwirda's inpoly
%i(i) = utils.inpoly([real(p(i(:))),imag(p(i(:)))], [real(v),imag(v)]);

%   MEX interface to Wm. Randolph Franklin's C code:
i(i) = logical(utils.inpolyc(p(i), v));   % converts int to logical
