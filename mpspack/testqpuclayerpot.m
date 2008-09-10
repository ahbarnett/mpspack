% test photonic bands QP unit cell layer pot sticking-out scheme
% barnett 8/6/08, convergence 8/10/08, checked for new evalcopies 9/8/08

clear classes
test = 'c';           % b = basis objects, c = convergnece
verb = 2;
k = 10;
uc = qpunitcell(1, 0.5+1i, k, 20);
uc.setbloch(-1, -1);               % furthest point from origin in Brillouin
if test=='b'     % ................. basis objects basic test
  if verb>1          % simple basis eval
    b = qpuclayerpot(uc, 'B', 's', k);
    uc.plot; b.showgeom; axis equal;
    sig = ones(b.Nf, 1);              % toy density
    dx = 0.04; gx = -2:dx:2; gy = -1:dx:1;
    [xx yy] = meshgrid(gx, gy); zz = xx(:)+1i*yy(:); A = b.eval(pointset(zz));
    u = reshape(A*sig, [numel(gy) numel(gx)]);
    showfield(gx, gy, u, [], 'qpuclayerpot.eval SLP test');
  end
  % time to see if reusing data is correct...
  tic; [A d] = b.evalunitcellcopies(pointset(zz)); fprintf('eval in %.2g s\n', toc)
  opts.data = d;
  tic; [A d] = b.evalunitcellcopies(pointset(zz), [], opts); fprintf('eval (pre-stored) %.2g s\n', toc)
  if 0                     % speed test for reusing data
    profile clear; profile on;
    for j=1:100, [A d] = b.evalunitcellcopies(pointset(zz), [], opts); end
    profile off; profile viewer
  end
  uc.addqpuclayerpots;
  if verb>1, A = uc.evalbases(pointset(zz));
    sig = ones(4*uc.bas{1}.Nf, 1);
    u = reshape(A*sig, [numel(gy) numel(gx)]);
    showfield(gx, gy, u, [], 'full addqpuclayerpot uc evalbases');
  end
  tic; Q = uc.evalbasesdiscrep; fprintf('eval discrep in %.2g s\n', toc)
  uc.setbloch(0);
  tic; Q = uc.evalbasesdiscrep; fprintf('eval discrep (pre-stored) %.2g s\n', toc) % should be much faster (using stored)
  
elseif test=='c'       % .............. convergence test

  b = mfsbasis(1.0+0.1i, [], [], k);     % make a source as a single far charge
  uc.addqpuclayerpots(k);
  pr = bvp(uc);                       % set up dummy problem w/ interior grid
  o.dx = .02; o = pr.gridboundingbox;
  gx = o.bb(1):o.dx:o.bb(2); gy = o.bb(3):o.dx:o.bb(4);
  [xx yy] = meshgrid(gx, gy); zz = xx + 1i*yy; ii = find(uc.inside(2*zz));
  p = pointset(zz(ii)); us = NaN*zz; us(ii) = b.eval(p);  % eval src field grid
  
  Ns = 5:5:40;
  rs = 0*Ns; idns = 0*Ns; tfs = 0*Ns; ts = 0 *Ns;     % resid, int norm, times
  for i=1:numel(Ns) % ====== loop over N
    fprintf('N = %d:\n', Ns(i))
    uc.seg.requadrature(Ns(i)); uc.setupbasisdofs;  % N quadr pts per UC wall
    tic; Q = uc.evalbasesdiscrep; tfs(i) = toc;
    fprintf('\tfill Q in %.2g s\n', tfs(i))
    fprintf('\tcond(Q) = %.3g\n', cond(Q))
    di = b.evalunitcelldiscrep(uc);        % eval discrep col vec
    tic; co = -Q\di; ts(i) = toc;
    rs(i) = norm(Q*co + di);
    fprintf('\tcoeff nrm = %.2g, resid l2 error = %.2g\n',norm(co),rs(i))
    pr.co = co; pr.setupbasisdofs;         % pass in co to dummy bvp
    [u di] = pr.pointsolution(p);          % only eval u at interior pts p
    idns(i) = norm(u+us(ii))*o.dx;
    fprintf('\tL2 int domain norm (u+us) = %.2g\n', idns(i))
  end             % =========
  figure; semilogy(Ns, [rs; idns; tfs; ts], '+-');
  title(sprintf('convergence, QP sticking-out layers, L2 err norm in central 1/4 of unit cell, k=%g', k));
  legend('l2 resid', 'L2 int domain u+us', 't_{fill} (s)', 't_{solve} (s)');
  xlabel('N per uc side'); axis([min(Ns) max(Ns) 1e-16 1])
  %print -depsc2 conv_qpuclayerpot.eps
end
