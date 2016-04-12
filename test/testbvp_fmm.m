% test MPSpack iterative Helmholtz BVP demo with FMMLIB2D/LP2D
% Barnett 3/21/11, tweaked 10/17/12
clear; % see set-up notes in Manual, for FMM and for LP2D codes.

N=1e3; k=10; verb=1;   % interior Dirichlet BVP for Helmholtz; try N=1e4...
s = segment.smoothstar(N,0.3,7); % actually has tight inside curvature, tricky!
d = domain(s, 1);   % interior
d.addlayerpot(s, [0 1], struct('quad','a','ord',8)); % DLP, Alpert quad
f = @(z) besselh(0,k*abs(z-2-3i));  % solution function giving bdry data
s.setbc(-1, 'd', [], @(t) f(s.Z(t)));
p = bvp(d); p.setoverallwavenumber(k);

%o.FMM = 0; o.meth='direct'; % 0 for dense solve (1 requires LP2D Alpert codes)
o.FMM = 1; o.meth = 'iter'; % FMM w/ GMRES for iterative soln, for large N
fprintf('testing N=%d; please wait about %g min...\n', N, N/60000); 
tic; p.solvecoeffs(o); fprintf('solve done in %.3g sec\n', toc)

% select o.FMM 0 or 1 here to use FMM for potential evaluation:
figure; o.FMM=1; o.dx=0.01; o.comparefunc=f; o.logabs=1; tic; p.showsolution(o);
fprintf('solution difference field eval in %.3g sec\n', toc)
% even at N=1e3, dx=0.01, FMM is 20x faster than direct evaluation!
