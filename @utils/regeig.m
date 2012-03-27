function [d V] = regeig(F, G, opts)
% REGEIG - solve Fix-Heiberger (Vergini) style regularized generalized EVP
%
% [d V] = REGEIG(F, G)
%
% [d V] = REGEIG(F, G, opts)

if nargin<3, opts = []; end
if ~isfield(opts, 'eps'), opts.eps = 1e-14; end   % defaults
if ~isfield(opts, 'v'), opts.v = 0; end

N = size(F,1);
[V,D] = eig(G); D = diag(D); i = find(D>opts.eps*max(D)); % indices to keep
clear G
r = numel(i); % rank
if opts.v, fprintf('\tN=%d, eps = %g, rank = %d\n', N, opts.eps, r); end
V = V(:,i) .* repmat(1./sqrt(D(i)).', [N 1]);  % compute V*L (rescale evecs)
tF = V'*F*V;               % projects out numerical Nul(G)
clear F
tF = (tF + tF')/2;         % make Hermitian
if nargout>1               % want eigvecs too
  [W, L] = eig(tF);
  d = diag(L);
  V = V*W;                 % transform evecs back to gevecs of (F,G)
  %figure; imagesc(log10(abs(V'*G*V - eye(r)))); title('V''GV');colorbar;caxis ([-16 0]);
  %figure; imagesc(log10(abs(V'*F*V - diag(l)))); title('V''FV');colorbar;caxis([-16 0]);
else                       % just eigvals
  d = eig(tF);
end
