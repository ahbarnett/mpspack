function Fh = spectralfiltermatrix(s,k,opts)
% SPECTRALFILTERMATRIX  fill M*M boundary matrix applying func of surf Laplacian
%
% Fh = spectralfiltermatrix(s,k,opts) returns square matrix spectrally-accurate
%   discretization of a certain function of the surface Laplacian Delta, ie
%
%     F_h = max [ (1 - Delta)_+^{1/2}, c/k^{1/3} ]
%
%   where (...)_+ indicates the positive part, c is a build-in constant, and
%   k is the overall wavenumber. This is a helper for filtered DtN method
%   for Neumann MPS bounds. See our paper arxiv:1512.04165
%
% Inputs:
%  s - segment
%  k - wavenumber
%  opts.verb - verbosity
%
% Switched to method that uses multiple of identity to handle the high-freq
% limit, and use Fourier projection for only the difference from this limit.

% Barnett 12/8/15. Brought in as @evp method 4/12/16.

if nargin==0, test_spectralfiltermatrix; return; end
if nargin<3, opts = []; end
if ~isfield(opts,'verb'), opts.verb = 0; end

L = sum(s.w);  % perim
M = numel(s.x);
arcl = real(ifft(-1i*fft(s.speed).*[0 1./(1:M/2-1) 0 -1./(M/2-1:-1:1)]'));

% spectral approx to sampled cumulative arclength
arcl = arcl'/2/pi + L*s.t'; % add back in growing component (const irrelevant)
maxm = ceil(M/4);  % can't make this big since get osc crap
%maxm = ceil(M/2)-1;
ns = -maxm:maxm; P = exp((-2i*pi/L)*ns'*arcl); % ns=freqs. Dense DFT mat
Pt = P'; PW = P.*repmat(s.w/L,[numel(ns) 1]); % Pt takes Fou coeffs to vals on pO
%norm(Pt*PW)  % close to Id, so close-ish to 1.
if opts.verb
  fprintf('test PW.1 should be zero apart from mode m=0 which is 1:\n')
  sum(PW,2)
end

% now stuff specific to filter func...
xin = 2*pi/L/k*ns;  % freq xi grid (wavenumber scaled so k=1)
c = 1.0;
if 0        % crude way that kills freqs > nmax (bad for use in tension)
  gk = max(real(sqrt(1-xin.^2)),c*k^(-1/3));   % hard h^{1/3} cutoff
  fk = 1./gk;
  Fh = Pt*(diag(fk)*PW);   % discrete inverse op, ie F_h
else        % use Id as the h^{1/3} shift, so high freqs handled correctly
  cut = c*k^(-1/3);   % cutoff and use as shift
  fk = 1./max(real(sqrt(1-xin.^2)),cut) - 1/cut;   % shift
  Fh = (1/cut)*eye(M) + Pt*(diag(fk)*PW);
end
Fh = real(Fh);         % matrix supposed to be real-valued
%%%%%

function test_spectralfiltermatrix
M=400;
s = segment.smoothnonsym(M, 0.3, 0.2, 3);
k = 30;
f = sin(k*real(s.x));  % some osc func going up to freqs k, but far from Nyq
Fh = evp.spectralfiltermatrix(s,k);
g = Fh*f;
figure; plot(cumsum(s.w), [f g], '.-');
%keyboard
