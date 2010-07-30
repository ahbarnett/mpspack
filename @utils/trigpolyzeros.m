function [r derr] = trigpolyzeros(F, opts)
% TRIGPOLYZEROS - return list of zeros of complex 2pi-periodic trig poly
%
% r = trigpolyzeros(F) returns list r of roots in (-pi, pi] of 2pi-periodic
%   function given by trigonometric polynomial with Fourier series coeff vector
%   F. The ordering of F is as returned by fftshift(fft(fftshift(f))), where f
%   samples the 2pi-periodic function at -pi+2*pi*(0:N-1)/N, ie complex
%   exponentials from frequency -N/2 to N/2-1, where N=length(F). The trig
%   poly's highest freq is chosen to be real (cos(N/2 t)) as in Trefethen,
%   Spectral Methods in Matlab book. Effort scales as O(N^3)
%
% r = trigpolyzeros(F, opts) passes in options including the following:
%   opts.tol: tolerance for how small in abs value must be to count as a zero
%             (default 1e-8)
%   opts.real: if true, assumes coeffs come from real func, keeps only UHP roots
%             (default false)
%
% [r derr] = trigpolyzeros(...) also returns distances of roots from unit circle

% (C) 2009, Alex Barnett.
  if nargin<2, opts = []; end
  if ~isfield(opts, 'tol'), opts.tol = 1e-8; end
  if ~isfield(opts, 'real'), opts.real = 0; end
  F = F(:);                                  % make col vector
  r = roots([F(1)/2; F(end:-1:2); F(1)/2]);  % degree-doubling a la Boyd 2002
  %figure; plot(r, '+'); hold on; plot(exp(1i*(0:.01:2*pi)),'-r'); axis equal
  derr = abs(abs(r)-1);
  % keep only close to unit circle, if real only those in upper half plane...
  ii = find(derr<=opts.tol & (~opts.real | imag(r)>=0));
  r = angle(r(ii));
  derr = derr(ii);

