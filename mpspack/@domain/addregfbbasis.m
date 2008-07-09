% ADDREGFBBASIS - create a regular Fourier-Bessel basis object in a domain
%
%  d = ADDREGFBBASIS(d, origin, N, k, opts) creates a regular FB basis
%   object within a domain object d, returning an updated copy of d.
%   Other arguments are as in REGFBBASIS
%
% See also: REGFBBASIS

function d = addregfbbasis(d, varargin)

d.bas  = {d.bas{:}, regfbbasis(varargin{:})}; % append cell arr of basis handles

if numel(varargin)>3
  d.k = varargin{4};                        % resets domain wavenumber
end
