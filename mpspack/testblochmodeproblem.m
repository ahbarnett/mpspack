% test Bloch mode problem class
% barnett 8/14/08
% test = 'c': convergence tests:
%             'cd' use fixed density for obstacle and QP LP basis, give some u
%             'cp' fix density on obstacle, then solve for QP field, give one u
%             'cs' solve for min sing val of whole matrix
% test = 's': sweep test, sweeps smallest sing vals over alpha.

clear classes
test = 's';
verb = 0;
k = 10; ns = 100;
s = scale(segment.smoothstar(50, 0.3, 3), 0.2);
e = domain([], [], s, -1);
e.addlayerpotbasis(s, 'd', k);
s.setbc(1, 'd');
uc = qpunitcell(1, 0.5+ 1i, k, 20);   % M=20 for sweep test ('g' def quadr)
uc.addqpuclayerpots;      % choose QP basis set
o.sourcenei = 1; pr = blochmodeproblem(uc, e, [], o); % choose # nei (src,discr)
if verb, figure; pr.showbdry; end

if test(1)=='c' % .......................... convergence tests....
  uc.setbloch(exp(1i*pi*0.8264), exp(1i*pi*0));  % pi*0.8264 is close to mode
  Ns = 10:5:40;
  ss = NaN*zeros(ns, numel(Ns));
  for i=1:numel(Ns) % ====== loop over N
    uc.seg.requadrature(Ns(i));       % set N quadr points per unit cell wall
    %uc.clearbases; o.rescale_rad = 0.8; uc.addregfbbasis(0, Ns(i), k, o);
    %uc.addrpwbasis(Ns(i), k);        % drop in these lines to try other QP bas
    if test(2)=='d'        % use given density func
      M = pr.wholematrix;
      pr.wco = ones(size(M,2), 1); % density = 1 on obstacle and QP LPs
      u = pr.pointsolution(pointset([0.4 + 0.3i]));
      fprintf('N=%d: \tu(pO)=%.15g \tdiscrep=%.15g \tu(x)=%.15g\n', Ns(i), M(1,:)*pr.wco, M(end,:)*pr.wco, u)
    elseif test(2)=='p'                % solve for QP field (periodize obst src)
      pr.setupbasisdofs; pr.co = ones(pr.N, 1);
      tic; di = pr.fillobstodiscrep * pr.co; toc
      tic; Q = uc.evalbasesdiscrep; toc
      qpco = -Q\di;
      uc.setupbasisdofs; pr.wco = [pr.co; qpco];      % stack whole coeff vec
      u = pr.pointsolution(pointset([0.4 + 0.3i]));
      fprintf('\t\t\t\tu(x) = %.15g\n', u);
    elseif test(2)=='s'                % compute min sing val soln
      o = []; o.twonorm = 1; M = pr.wholematrix(o);
      [U S V] = svd(M); sing = diag(S); % sing = svd(M);
      n = min(ns,numel(sing)); ss(1:n,i) = sing(end-n+1:end); % save data
      fprintf('N = %d: min sing val = %.15g\n', Ns(i), min(ss(:,i)))
      pr.wco = V(:,end); pr.wholecoeffvectwonormreweight;
      u = pr.pointsolution(pointset([0.4 + 0.3i; 0.4 + 0.2i]));
      fprintf('\t\t\t\tu(x1)/u(x2) = %.15g\n', u(1)/u(2));
    end
  end               % ======
if verb, figure; imagesc(log10(abs(ss-1))); colorbar; title('log10 deviations of \sigma from 1'); end

elseif test=='s' %....................... sweep test
be = exp(1i*pi*0); phs = pi*(0.8234:0.001:0.8294);   % sweep phases
ss = NaN*zeros(ns, numel(phs));         % singular vals data (each col)
%profile clear; profile on;
for i=1:numel(phs)
  al = exp(1i*phs(i)); uc.setbloch(al, be); fprintf('alpha ang = %.4g pi\n',phs(i)/pi)
  o = []; o.twonorm = 1; M = pr.wholematrix(o);
  [U S V] = svd(M); sing = diag(S);
  n = min(ns,numel(sing)); ss(1:n,i) = sing(end-n+1:end);
  if 0 & abs(phs(i)-0.827*pi)<1e-3
    pr.wco = V(:,end); pr.wholecoeffvectwonormreweight;     % right sing vec
    opts.dx = 0.01; tic; [uN gx gy di] = pr.gridsolution(opts); toc;
    figure; imagesc(gx, gy, real(uN)); title('Bloch mode');
    c=caxis; caxis(max(c)*0.3*[-1 1]); axis equal tight; colorbar;
    set(gca,'ydir','normal');drawnow;
    [xx yy] = meshgrid(gx, gy); zz = xx+1i*yy;
    bw = exp(-1i*real(conj(uc.kbloch).*zz));    % conj of bloch wave part
    figure; imagesc(gx, gy, real(uN.*bw)); title('Periodic part of mode');
    c=caxis; caxis(max(c)*0.3*[-1 1]); axis equal tight; colorbar;
    set(gca,'ydir','normal');drawnow;
  end
end
%profile off; profile viewer
figure; plot(phs, ss', '-+'); xlabel('ang'); ylabel('\sigma_j(M)');
v=axis; v(3:4) = [0 1e-3]; axis(v);
title('min sing val, \alpha sweep (\beta=1, k=10)');
end

if verb, figure; imagesc(log10(abs(M))); colorbar; title('log_{10} |M|'); end
if verb, figure; imagesc(real(M)); colorbar; title('Re M'); end
if verb, figure; plot([real(diag(M)) imag(diag(M))]);
  title('diag M (Re & Im)'); end
%[U,S,V] = svd(M);

