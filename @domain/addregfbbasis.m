% ADDREGFBBASIS - create a regular Fourier-Bessel basis object in a domain
%
%  ADDREGFBBASIS(d, origin, N, opts) creates a regular FB basis
%   object within a domain object whose handle is d.
%   Other arguments are as in REGFBBASIS
%
% See also: REGFBBASIS

function addregfbbasis(d, varargin)

d.bas  = {d.bas{:}, regfbbasis(varargin{:})}; % append cell arr of basis handles

d.bas{end}.doms = d;                    % tell this basis it affects this domain
