% Test the Jfilter close-eval of QP LP scheme, to debug filling B matrix in
% testblochbyhand with Jfilter. Adapted from testblochbyhand.m
% Barnett 1/19/09

clear all classes
k = 10; Mu = 50; Ms = 40;
uc = qpunitcell(1, 0.5+1i, k, Mu); uc.buffer = 0;
uc.addqpuclayerpots;   % use this rather than following line for full test
%uc.addlayerpotbasis(uc.L, 'd', k); uc.addlayerpotbasis(uc.B, 'd', k); % simpler
uc.setupbasisdofs;
s = scale(segment.smoothstar(Ms, 0.3, 3), 0.2);
e = domain([], [], s, -1); o.dom = e;
b = e.addlayerpotbasis(s, [1i*k 1], k); b = b{1};

if 1
  tic; Q = uc.evalbasesdiscrep; toc
  o.nei = 1; tic; C = b.evalunitcelldiscrep(uc, o); toc
  nu = ones(size(s.x));                         % given density on segment
  qpdens = Q\(C*nu);                  % solve for qp lp densities vector
  %qpdens = 0*qpdens; qpdens(10) = 1; %qpdens(2*Mu+1:end) = 0;
else
  qpdens = ones(uc.N, 1);
end

for meth=1:2         % --------- loop over tests: 1=direct eval, 2=Jfilter eval
  bo = []; if meth==2, bo.Jfilter.M = 25; bo.Jfilter.rescale_rad = 0.5;
  end

  g = -2:0.02:2;
  [xx yy]=meshgrid(g,g); p=pointset(xx(:)+1i*yy(:)); disp('eval QP field...')
  tic; B = uc.evalbases(p, bo); toc;
  u = reshape(B * qpdens, size(xx));
  figure; imagesc(g, g, real(u)); set(gca, 'ydir','normal'); axis equal;
  uc.plot; xlabel x; ylabel y; title('QP field');
  if meth==1, unaive=u; c = max(abs(caxis)); end; caxis(c*[-1 1]);
end

figure; imagesc(g, g, log10(abs(u-unaive))); caxis([-16 0]); colorbar;
set(gca, 'ydir','normal'); axis equal; xlabel x; ylabel y;
title('QP field err log10');
