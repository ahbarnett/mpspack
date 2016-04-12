function [d V] = filteredDtNspectrum(p, k, o)
% FILTEREDDTNSPECTRUM - eigenvalues/eigenfunctions of filtered DtN bdry op
%
% d = filteredDtNspectrum(p, kstar) returns all eigenvalues d of the weighted
%  Dirichlet-to-Neumann operator Theta(kstar), at wavenumber kstar, used
%  for Neumann scaling-type method.
%
% Currently the k used for spectral filter can only be fixed as a hack.
%
% See NtDspectrum for other notes; some code from neuscaling.m

% Copyright (C) 2014, Alex Barnett, based on NtDspectrum

if numel(p.segs)~=1, error('evp object must contain exactly 1 segment!'); end
s = p.segs(1);            % get the one segment
M = numel(s.w); L = sum(s.w); % # pts, perimeter
arcl = real(ifft(-1i*fft(s.speed).*[0 1./(1:M/2-1) 0 -1./(M/2-1:-1:1)]'));
% spectral approx to sampled cumulative arclength
arcl = arcl'/2/pi + L*s.t'; % add back in growing component (const irrelevant)

if nargin<3, o = []; end  % process options - from NtDspectrum
if ~isfield(o, 'quad'), o.quad = 'm'; end   % default layerpot quadrature corrn
if ~isfield(o, 'cayley'), o.cayley = 0; end
eta = k;                  % inverse scale param in Cayley xform
if isfield(o, 'wei'), w = o.wei; else w = 1./real(conj(s.x).*s.nx); end % 1/x.n

N = numel(s.x);
HpD = eye(N)/2 + layerpot.D(k, s, [], o);               % 1/2 + D
S = layerpot.S(k, s, [], o);                            % S

if ~o.cayley
  DtN = inv(S) * HpD;
else, error('cayley not implemented'); end
clear Sw HpD

% set up Fourier filter on bdry; see neuscaling...
persistent PW Pt kicknull  % only compute once (indep of k)
ns = -N/4:N/4; %ns = -N/2+1:N/2;  % don't go out to all freqs where quadr inacc
if isempty(PW)
P = exp((-2i*pi/L)*ns'*arcl); % ns=freqs. Dense DFT mat
Pt = P'; PW = P.*repmat(s.w/L,[numel(ns) 1]); % Pt takes F coeffs to vals on pO
[U S V] = svd(PW); nulP = V(:,numel(ns)+1:end); % null PW
kicknull = 1e1*nulP*nulP'; clear U S V nulP % kick nullspace eigvals up to big
end

kfilt = k; %301;    % or k,  filter wavenumber
xin = 2*pi/L/kfilt*ns;  % freq xi grid (wavenumber scaled so k=1)
%Fk = max(real(sqrt(1-xin.^2)),0.5*kfilt^(-1/3)); % freq filter vs xi, default
Fk = max(real(sqrt(1-xin.^2)),0.1*kfilt^(-1/3)); % freq filter vs xi
%Fk = Fk .* (1 - (1-1e-2)*(abs(ns)>N/4));     % truncate in freq, so iF large
iF = Pt*(diag(Fk.^-1)*PW) + kicknull;  % high freqs kicked to big EVs
%iFk = max(real(1./sqrt(1-xin.^2))); iFk(isinf(iFk)) = 0; % inv of freq filter
%iF = Pt*(diag(iFk)*PW);  % fails?
wDtN = iF * (DtN * (iF .* repmat(w(:).', [N 1])));  % A^-1 * DtN * A^-1 * xn^-1

if nargout==1
  d = eig(wDtN);     % eigenvalues only: dense
else
  [V D] = eig(wDtN); % eigenvectors also: dense
  d = diag(D);
end
