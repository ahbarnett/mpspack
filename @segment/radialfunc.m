function s = radialfunc(M, fs, varargin)
% RADIALFUNC - closed segment from radial function r=f(theta) and derivatives
%
%  s = radialfunc(M, {f fp}) generates a smooth closed radial function segment
%   with M discretization pts, using periodic trapezoid quadrature. f and fp
%   are function handles to the radius function f(t) and its derivative, for
%   angle t in [0,2pi].
%
%  s = radialfunc(M, {f fp fpp}) also includes 2nd-derivative, enabling
%   curvature information to be created in the segment (for spectral quadr)
%
% Copyright (C) 2008, 2009, Alex Barnett, Timo Betcke

f = fs{1};
Z = @(s) exp(2i*pi*s).*f(2*pi*s);   % note conversion from 0<s<1 to 0<t<2pi
fp = fs{2};
if numel(fs)<=2
  s = segment(M, {Z, @(s) 2*pi*(1i*Z(s) + exp(2i*pi*s).*fp(2*pi*s))}, 'p');
else
  Zp = @(s) 2*pi*(1i*Z(s) + exp(2i*pi*s).*fp(2*pi*s));
  fpp = fs{3};
  Zpp = @(s) 4*pi^2*((fpp(2*pi*s) + 2i*fp(2*pi*s)).*exp(2i*pi*s) - Z(s));
  s = segment(M, {Z, Zp, Zpp}, 'p', varargin{:});
end
% Notice that all this routine did was manipulate function handles, no evals.
