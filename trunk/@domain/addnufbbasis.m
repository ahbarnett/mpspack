function addnufbbasis(d, varargin)
% ADDNUFBBASIS - create a fractional-order Fourier-Bessel basis set (corner exp)
%
%   ADDNUFBBASIS(d, origin, nu, offset, branch, N, opts) creates a basis of fractional-
%   order Fourier-Bessel functions appropriate for expansion of the Helmholtz
%   equation in a wedge of angle pi/nu > 0. The orders are nu*(1:N) (for sine
%   angular functions) or nu*(0:N) (for cosine angular functions). For a
%   full description see NUFBBASIS
%
%  See also: NUFBBASIS

% Copyright (C) 2008, 2009, Alex Barnett, Timo Betcke

d.bas  = {d.bas{:}, nufbbasis(varargin{:})}; % append cell arr of basis handles

d.bas{end}.doms = d;                    % tell this basis it affects this domain
