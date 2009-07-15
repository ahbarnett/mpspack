% try basic interior BVP 'by hand' in order to see best way to do bvp class.
% also see testbasisdomain.m
% barnett 7/10/08. Now somewhat obsolete, since bases no longer inside domains.

clear all classes
verb = 1;
M = 40; s = segment.polyseglist(M, [1 1i exp(4i*pi/3)]);  % CCW tri from INCL
d = domain(s, 1);
k = 10; N = 30;
opts.real = 1; d.addrpwbasis(N, k, opts);

% set up BCs...
%dirfunc = @(x) besselj(0, k*abs(x - 0.3 + .2i));   % boundary Dirichlet data
dirfunc = @(x) bessely(0, k*abs(x - 1.3 + 1.2i));   % boundary Dirichlet data
for j=1:numel(s)
  s(j).a = [0 1]; s(j).f = dirfunc(s(j).x);     % set BC for each seg
end
f = vertcat(s.f); if verb>1, figure; plot([real(f) imag(f)]); end % view data
% NB, simple since only M discrep vals per seg (no continuity conditions)

% fill A matrix for each domain (here there's only 1)
A = [];
for j=1:numel(s)
  seg = s(j);
  As = d.evalbases(seg);
  A = [A; As];
end

alpha = A \ f;   % solve least-squares

% show and check solution on interior grid...
dx = 0.01; [zz ii g h] = d.grid(dx);
fd = dirfunc(zz);                      % true solution
Ad = d.evalbases(pointset(zz));        % evaluation matrix for this domain grid
uNd = Ad * alpha;                      % our solution eval on grid
r = dx * norm(uNd - fd);               % L2 interior error norm (est on grid)
fprintf('coeff l2 norm = %g, residual L2 norm in domain = %g\n', norm(alpha),r);
if verb
  uN = NaN*zeros(size(ii)); uN(ii) = uNd;         % write uN to a grid
  figure; subplot(1,2,1); imagesc(g, h, uN); colorbar; axis equal tight;
  title('uN');
  e = NaN*zeros(size(ii)); e(ii) = uNd - fd;      % write pointwise err to grid
  subplot(1,2,2); imagesc(g, h, e); colorbar; axis equal tight;
  title('uN - u_{exact}');
end

% convergence plot ? etc...
