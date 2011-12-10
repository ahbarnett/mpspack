function [d V] = NtDspectrum(p, k, o)
% NTDSPECTRUM - return eigenvalues/eigenfunctions of (weighted) NtD bdry op
%
% d = NtDspectrum(p, kstar) returns all eigenvalues d of the weighted
%  Neumann-to-Dirichlet operator Theta(kstar), at wavenumber kstar.
%
% [d V] = NtDspectrum(p, kstar) also returns all eigenvectors
%
% [d ...] = NtDspectrum(p, kstar, opts) controls options including:
%   opts.quad: quadrature correction scheme for segment later-potentials
%              (default is 'm', Kress scheme; other option is 'a', Alpert)
%   opts.ord: order for Alpert quadrature correction (4,8,16, etc).
%   opts.cayley: if true (default), use Cayley transform, otherwise naive way.
%   opts.wei: if present, overrides the vector of weights 1/(x.n) at nodes
%
% Notes:
% 1) Currently only support a domain bounded by a single closed segment
% 2) Dense linear algebra is used.
% 3) Typical imag part of returned eigenvalues serves as a discretization error
%    overall estimate.
%
% See also: EVP

% Copyright (C) 2011, Alex Barnett

if numel(p.segs)~=1, error('evp object must contain exactly 1 segment!'); end
s = p.segs(1);            % get the one segment

if nargin<3, o = []; end  % process options
if ~isfield(o, 'quad'), o.quad = 'm'; end   % default layerpot quadrature corrn
if ~isfield(o, 'cayley'), o.cayley = 1; end
eta = k;                  % inverse scale param in Cayley xform
if isfield(o, 'wei'), w = o.wei; else w = 1./real(conj(s.x).*s.nx); end % 1/x.n

N = numel(s.x);
HpD = eye(N)/2 + layerpot.D(k, s, [], o);               % 1/2 + D
Sw = layerpot.S(k, s, [], o) .* repmat(w(:).', [N 1]);  % S (x.n)^{-1}

if ~o.cayley
  wNtD = inv(HpD) * Sw;
else
  wNtD = inv(HpD - 1i*eta*Sw) * (-HpD - 1i*eta*Sw);
end
clear Sw HpD

if nargout==1
  d = eig(wNtD);     % eigenvalues only: dense
else
  [V D] = eig(wNtD); % eigenvectors also: dense
  d = diag(D);
end

if o.cayley, d = (1i/eta) .* (1+d)./(1-d); end  % undo Cayley xform each eigval
