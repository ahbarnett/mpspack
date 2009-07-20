function addregfbbasis(d, varargin)
% ADDREGFBBASIS - create a regular Fourier-Bessel basis object in a domain
%
%  ADDREGFBBASIS(d, origin, N, opts) creates a regular FB basis
%   object within a domain object whose handle is d.
%   The rest of the argument list is discussed in REGFBBASIS
%
% See also: REGFBBASIS
%
% Copyright (C) 2008, 2009, Timo Betcke, Alex Barnett


d.bas  = {d.bas{:}, regfbbasis(varargin{:})}; % append cell arr of basis handles

d.bas{end}.doms = d;                    % tell this basis it affects this domain
