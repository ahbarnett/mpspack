% ADDRPWBASIS - create a real plane wave basis object in a domain
%
%  ADDRPWBASIS(d, N, opts) creates a real-valued plane wave basis
%   object within a domain object whose handle is d.
%   Other arguments are as in RPWBASIS
%
% See also: RPWBASIS

function addrpwbasis(d, varargin)

d.bas  = {d.bas{:}, rpwbasis(varargin{:})}; % append cell arr of basis handles
