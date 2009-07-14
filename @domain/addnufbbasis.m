% ADDNUFBBASIS - create a regular Fourier-Bessel basis object in a domain
%
%  ADDNUFBBASIS(origin, nu, offset, branch, N, opts) creates an irregular FB
%   basis object within a domain object whose handle is d.
%
% See also: NUFBBASIS

function addnufbbasis(d, varargin)

d.bas  = {d.bas{:}, nufbbasis(varargin{:})}; % append cell arr of basis handles
