function f = smoothedstep(t, opts)
% f = smoothedstep(t) returns a set of values of a smoothed step function
%  (roll-off, or, really roll-up) given array of arguments t. The returned f
%  has the same size as t.
%  The default function is integral of a Gaussian: f(t) = (1+erf(chi(t-1/2)))/2
%  The graph has inversion symmetry about the point (1/2,1/2).
%
% f = smoothedstep(t, opts) allows options to be set, including:
%  opts.chi : determines steepness parameter chi (default 10)

% Copyright (C) 2011, Alex Barnett
if nargin<2, opts = []; end
if ~isfield(opts, 'chi'), opts.chi = 10.0; end
f = (erf(opts.chi * (t-1/2)) + 1)/2;
