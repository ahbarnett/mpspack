% ADDREGFBBASIS - create a regular Fourier-Bessel basis object in a domain
%
%  ADDREGFBBASIS(d, origin, N, k, opts) creates a regular FB basis
%   object within a domain object whose handle is d.
%   Other arguments are as in REGFBBASIS
%
% See also: REGFBBASIS

function addregfbbasis(d, varargin)

d.bas  = {d.bas{:}, regfbbasis(varargin{:})}; % append cell arr of basis handles

if numel(varargin)>2
  d.k = varargin{3};                        % resets domain wavenumber
end
