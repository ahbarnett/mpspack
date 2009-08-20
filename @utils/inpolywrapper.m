function i = inpolywrapper(p, v)
% INPOLYWRAPPER - logical array of if numbers p are in polygon v (various codes)
%
% i = INPOLYWRAPPER(p, v) where p and v are column vectors of complex numbers,
%   returns logical array i of size same as p, giving whether each point p is
%   within simply-connected polygon with vertices v.
%
% Chooses between one of 3 methods:
%   * Matlab's native inpolygon
%        inpolygon_vec is very slow, 0.5us per point per vertex
%   * Darren Engwirda's inpoly
%        claims to have some intelligent sorting of points in p, but I find
%        only 1.25 times speed of matlab
%   * MEX interface to Wm. Randolph Franklin's C code
%        brute forces it and loops over points in p, seems to be 70 times
%        faster than matlab, ie only 8 ns per point per vertex!
%
% See also: INPOLYGON, TEST/TESTINPOLYWRAPPER, UTILS/INPOLYC

% Copyright (C) 2008, 2009, Timo Betcke, Alex Barnett

% uncomment one of the following calls...

%   Matlab's native inpolygon
i = inpolygon(real(p), imag(p), real(v), imag(v));

%   Darren Engwirda's inpoly
%i = utils.inpoly([real(p(:)),imag(p(:))], [real(v),imag(v)]);

%   MEX interface to Wm. Randolph Franklin's C code:
% i = logical(utils.inpolyc(p, v));   % converts int to logical
