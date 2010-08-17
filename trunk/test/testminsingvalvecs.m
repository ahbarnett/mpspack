% test routine for iterative alg util.minsingvalvecs on square matrix svd
% barnett 8/16/10

N = 700; iscomplex = 1;                         % choose size, and 0 or 1 here
A = rand(N) - 0.5 + iscomplex*1i*(rand(N)-0.5);
fprintf('wait for dense order-%d svd...\n',N)
tic; [U S V] = svd(A); v0 = V(:,end); u0 = U(:,end); s0 = S(end,end); toc
if v0(1)<0, v0 = -v0; u0 = -u0; end   % ensure positive real first entry
fprintf('s true=%.16g\n',s0)

disp('testing minsingvalvecs...')
for i=1:4, o=[]; if i>2; o.tol=1e-6; fprintf('opts.tol=%.3g...\n',o.tol), end
  tic; [u s v info] = utils.minsingvalvecs(A,o); toc
  fprintf('s=%.16g (err=%.3g), info.flag=%d, info.its=%d\n', s, s-s0, ...
          info.flag, info.its)
  fprintf('norm(A*v-s*u)=%.3g, norm(A''*u-s*v)=%.3g\n',...
          norm(A*v-s*u), norm(A'*u-s*v))
  fprintf('u errnorm=%.3g, v errnorm=%.3g, v proj err=%.3g\n\n', norm(u-u0),...
          norm(v-v0), norm(v-(v0'*v)*v0)) % subspace angle error (ignores phase)
end
