% test Helmholtz-FMM'able basis set evaluations: MFS, layerpot. Barnett 3/11/11
clear; verb = 0;
dx = 0.01; g = -1:dx:1; d = domain(); d.k = 10;   % display grid, wavenumber
basis = 's';               % which basis type to test: m, s, or d
N = 100; co = randn(N,1);  % basis size and coefficients
if basis=='m'  % MFS basis
  type = 3; etas = [inf, 0, -d.k]; opts.eta = etas(type); opts.real = 0;
  a=0.5; b = mfsbasis({@(t)a*exp(2i*pi*t),@(t)2i*pi*a*exp(2i*pi*t)}, N, opts);
else           % layerpot basis
  s = scale(segment.smoothstar(N,0.3,5), 0.501); % curve, keep from grid
  b = layerpot(s, basis);  % 's' or [1 0] for single, 'd' or [0 1] for double
end
b.doms = d; d.bas{1} = b;  % associate basis with domain

[xx yy] = meshgrid(g, g); M = numel(xx);
p = pointset([xx(:) + 1i*yy(:)], ones(size(xx(:)))); % set up n for x-derivs
fprintf('evaluating basis...\n')
tic; [A Ax Ay] = b.eval(p); t=toc;
fprintf('\t%.3g direct evals (& x,y derivs) in %.3g secs = %.3g evals/sec\n',...
        numel(A), t, numel(A)/t)
u = A*co; ux = Ax*co; uy = Ay*co;
tic; [v vx vy] = b.evalFMM(co, p); t=toc;     % do FMM
fprintf('\tHelmholtz FMM done %.3g secs = %.3g points (src+targ) per sec\n',...
        t, (N+M)/t);
fprintf('l2-errors: u = %.3g, ux = %.3g, uy = %.3g (large if s-t close)\n', ...
        norm(v-u), norm(vx-ux), norm(vy-uy))   % could test s.dist(p) close

if verb, figure; subplot(2,3,1); w=real(reshape(u,size(xx)));
  imagesc(g, g, w); set(gca,'ydir','normal'); axis equal tight;
  utils.goodcaxis(w); title('direct summation Re[u]'); drawnow;
  subplot(2,3,2); w=real(reshape(ux,size(xx)));
  imagesc(g, g, w); set(gca,'ydir','normal'); axis equal tight;
  utils.goodcaxis(w); title('direct summation Re[u_x]'); drawnow;
  subplot(2,3,3); w=real(reshape(uy,size(xx)));
  imagesc(g, g, w); set(gca,'ydir','normal'); axis equal tight;
  utils.goodcaxis(w); title('direct summation Re[u_y]'); drawnow;
  subplot(2,3,4); z=reshape(v-u,size(xx)); imagesc(g, g, log10(abs(z)));
  caxis([-16 0]); axis equal tight; title('log_{10} error in u'); %colorbar;
  subplot(2,3,5); z=reshape(vx-ux,size(xx)); imagesc(g, g, log10(abs(z)));
  caxis([-16 0]); axis equal tight; title('log_{10} error in u_x'); %colorbar;
  subplot(2,3,6); z=reshape(vy-uy,size(xx)); imagesc(g, g, log10(abs(z)));
  caxis([-16 0]); axis equal tight; title('log_{10} error in u_y'); %colorbar;
  
  % test problem class plotting...
  pr = bvp(d); pr.setupbasisdofs; pr.co = co; o.bb = 3*[-1 1 -1 1];
  figure; o.FMM = 1; tic; [v] = pr.showsolution(o); toc
  figure; o.FMM = 0; tic; [u] = pr.showsolution(o); toc
  fprintf('l2-error in pointsolution = %.3g\n', norm(v(:)-u(:)))
  figure; semilogy(abs(v(:)-u(:))); title('errors at gridpoints')
end
  
if 0 & basis~='m'  % LP so test self-interaction...
  N = 1e3; s.requadrature(N); co = randn(N,1);  % make more beefy test
  b.quad = 'a'; b.ord = 16;       % choose quadrature scheme for layerpot
  d = domain(s, 1); b.doms = d; d.bas{1} = b; % choose interior domain (for JRs)
  d.k = 10; o.dom = d;  % wave#, pass in as opt to eval
  fprintf('evaluating self-interaction of layerpot basis...\n')
  tic; [A An] = b.eval(s, o); t=toc;
  fprintf('\t%.3g direct evals (& n-derivs) in %.3g secs = %.3g evals/sec\n',...
          numel(A), t, numel(A)/t)
  u = A*co; un = An*co;
  tic; [v vn] = b.evalFMM(co, s, o); t=toc;     % do FMM
  fprintf('\tHelmholtz FMM done %.3g secs = %.3g pts (src+targ) per sec\n',...
          t, (N+M)/t);
  fprintf('l2-errors: u = %.3g, un = %.3g\n', norm(v-u), norm(vn-un))
end
