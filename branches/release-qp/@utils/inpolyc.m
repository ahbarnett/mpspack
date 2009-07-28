% INPOLYC - MEX interface to fast inpolygon implementation code by W.R.Franklin
%
%  i = inpolyc(p, v) returns a column vector with entries 0 or 1, of int32 type,
%   specifying whether each point in the complex column vector p is inside the
%   polygon with vertices given in the complex column vector v.
%
%   p and v may be of real rather than complex type, in which case zero imag
%   part is assumed.
%
% See also: UTILS.INPOLYWRAPPER, TEST/TESTINPOLYWRAPPER

% Copyright (C) 2008, 2009, Timo Betcke, Alex Barnett

% Original C-code acknowledgments: (see inpolyc.c for full statement)
% Copyright (c) 1970-2003, Wm. Randolph Franklin
