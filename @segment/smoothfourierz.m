function s = smoothfourierz(N, z, nterms, tol)
% SMOOTHFOURIERZ - general complex Fourier series closed segment
%
%  s = SMOOTHFOURIERZ(N, z) generates a smooth closed segment defined by a
%   numel(z) term Fourier series, passing through all points z in the complex
%   plane. numel(z) must be even for now. N defines the number of points on the
%   new curve.
% 
%  s = SMOOTHFOURIERZ(N, z,nterms) limits the Fourier series to nterms terms.
%   If nterms<numel(z), then the curve does not pass through the points z
%   exactly, rather, approximates them. nterms even for now.
%
%  s = SMOOTHFOURIERZ(N, z,nterms,tol) uses a soft roll-off of Fourier
%   coeffs inducing at most tol error in the low-freq coeffs, without Gibbs
%   ringing.
%
% see: test/testsmoothfourierz

% Copyright (C) 2020 Alex Barnett

n = numel(z);
if mod(n,2), error('numel(z) must be even'); end
if nargin<3, nterms = n; end
if mod(nterms,2), error('nterms must be even'); end
zhat = fft(z(:))/n;
if nargin>3           % smooth rolloff
  wid = nterms/2 / sqrt(log(1/tol));
  rolloff = erfc((-nterms/2:nterms/2-1)/wid)/2;
  rolloff = rolloff(1:n/2);
  zhat = zhat .* [rolloff rolloff(end:-1:1)].';
  nterms = min(2*nterms, n);
end
zhat = [zhat(1:nterms/2); zhat(end-nterms/2+1:end)];   % keep lowest nterms
s = segment(N,{@(t) fourierZ(zhat,t), @(t) fourierZp(zhat,t), @(t) fourierZpp(zhat,t)},'p');

% analytic formulae for a Fourier segment --------------
function z = fourierZ(zhat,t)     % must work on vector of t's
t = 2*pi*t;
N = numel(zhat);  % even
z = 0*t;
for k=0:N/2
  z = z + zhat(k+1)*exp(1i*k*t);
end
for k=-N/2:-1
  z = z + zhat(k+1+N)*exp(1i*k*t);
end

function zp = fourierZp(zhat,t);  % deriv func Z'
N = numel(zhat);
zp = 2*pi*fourierZ(zhat.*[0 1i*(1:N/2-1) 0 1i*(-N/2+1:-1)].', t);

function zpp = fourierZpp(zhat,t);  % deriv func Z''
N = numel(zhat);
zpp = 2*pi*fourierZp(zhat.*[0 1i*(1:N/2-1) 0 1i*(-N/2+1:-1)].', t);
