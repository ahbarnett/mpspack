% test photonic bands QP unit cell layer pot sticking-out scheme
% barnett 8/6/08, convergence 8/10/08, checked for new evalcopies 9/8/08
% converted to new release-qp interface barnett 7/30/09

clear classes
test = 'c';           % b = basis objects, c = convergnece
verb = 1;    % 1 does tests, 2 shows figures of fields
bo = [];                % basis opts
%bo.Jfilter.M = 22; bo.Jfilter.rescale_rad = 0.5; % use J-filter for QPLP eval?
k = 10;
uc = qpunitcell(1, 0.5+1i, k, 20);
uc.setbloch(-1, -1);               % furthest point from origin in Brillouin

if test=='b'  % ........................... basis objects basic test
  dx = 0.04; gx = -2:dx:2; gy = -1:dx:1;
  [xx yy] = meshgrid(gx, gy); zz = xx(:)+1i*yy(:);
  b = qpuclayerpot(uc, 'B', 's');
  if verb>1          % simple basis eval
    uc.plot; b.showgeom; axis equal;
    sig = ones(b.Nf, 1);              % toy density
    A = b.eval(pointset(zz), bo); u = reshape(A*sig, [numel(gy) numel(gx)]);
    showfield(gx, gy, u, [], 'qpuclayerpot.eval SLP test');
  end
  % time to see if reusing data is correct...
  tic; [A d] = b.evalunitcellcopies(pointset(zz)); fprintf('eval in %.2g s\n', toc)
  opts.data = d;
  tic; [Ap d] = b.evalunitcellcopies(pointset(zz), [], opts); fprintf('eval (pre-stored) %.2g s\n', toc)
  fprintf('\terr prestored from cold-computed = %g\n', norm(A-Ap))
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
  %uc.setbloch(0);
  tic; Qp = uc.evalbasesdiscrep; fprintf('eval discrep (pre-stored) %.2g s\n', toc) % should be much faster (using stored)
  fprintf('\terr prestored from cold-computed = %g\n', norm(Q-Qp))    
  
elseif test=='c' %.......... convergence test: cancel u in UC via known discrep

  b = mfsbasis(pointset(1.0+0.1i));    % make a source as a single far charge
  b.doms = domain(); b.doms.k = k;                   % have it influence R^2
  uc.addqpuclayerpots;
  pr = bvp(uc);                       % set up dummy problem w/ interior grid
  o.dx = .02; o = pr.gridboundingbox;
  gx = o.bb(1):o.dx:o.bb(2); gy = o.bb(3):o.dx:o.bb(4);
  [xx yy] = meshgrid(gx, gy); zz = xx + 1i*yy;
  ii = find(uc.inside(2*zz));  % keeps only central 1/4 of UC shape
     % note with Jfilter may replace with (zz), ie whole interior (conv slow)
  p = pointset(zz(ii)); us = NaN*zz; us(ii) = b.eval(p);  % eval src field grid
  
  Ns = 5:5:40;
  rs = 0*Ns; idns = 0*Ns; tfs = 0*Ns; ts = 0*Ns;     % resid, int norm, times
  for i=1:numel(Ns) % ====== loop over N
    fprintf('N = %d:\n', Ns(i))
    uc.seg.requadrature(Ns(i)); uc.setupbasisdofs;  % N quadr pts per UC wall
    tic; Q = uc.evalbasesdiscrep; tfs(i) = toc; % don't use Jfilter to fill Q!
    fprintf('\tfill Q in %.2g s\n', tfs(i))
    fprintf('\tcond(Q) = %.3g\n', cond(Q))
    di = b.evalunitcelldiscrep(uc);        % eval discrep col vec
    tic; co = -Q\di; ts(i) = toc;
    rs(i) = norm(Q*co + di);
    fprintf('\tcoeff nrm = %.2g, resid l2 error = %.2g\n',norm(co),rs(i))
    u = uc.evalbases(p, bo) * co;   % unlike pr.pointsolution, lets bas opts in
    idns(i) = norm(u+us(ii))*o.dx;
    fprintf('\tL2 int domain norm (u+us) = %.2g\n', idns(i))
  end             % =========
  figure; semilogy(Ns, [rs; idns; tfs; ts], '+-');
  title(sprintf('convergence, QP sticking-out layers, L2 err norm in central 1/4 of unit cell, k=%g', k));
  legend('l2 resid', 'L2 int domain u+us', 't_{fill} (s)', 't_{solve} (s)');
  xlabel('N per uc side'); axis([min(Ns) max(Ns) 1e-16 1])
  %print -depsc2 conv_qpuclayerpot.eps
  % Note that conv rate (ex Jfilt) dep on how close the pointset p is to walls
end
