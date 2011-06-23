% test MPSpack iterative Helmholtz BVP with LP2D/HFMM2D
% Barnett 3/21/11
clear; addpath /home/alex/physics/leslie/gimbutas/lp2d/

% linsolve gmres iterapp crashes for FMM=1 ..... fix this!

N=5e2; k=10; verb=1;   % interior Dirichlet BVP for Helmholtz
s = segment.smoothstar(N,0.3,7); % actually has tight inside curvature, tricky!
d = domain(s, 1);   % interior
d.addlayerpot(s, [0 1], struct('quad','a','ord',8)); % DLP
f = @(z) besselh(0,k*abs(z-2-3i));  % solution function giving bdry data
s.setbc(-1, 'd', [], @(t) f(s.Z(t)));
p = bvp(d); p.setoverallwavenumber(k);
o.FMM = 1; o.meth = 'iter'; %o.meth = 'direct';
tic; p.solvecoeffs(o); fprintf('solve done in %.3g sec\n', toc)
figure; o.comparefunc=f; o.logabs=1; tic; p.showsolution(o);
fprintf('solution difference field eval in %.3g sec\n', toc)

