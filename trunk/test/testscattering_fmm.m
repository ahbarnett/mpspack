% test MPSpack iterative Helmholtz sound-soft scattering with LP2D/HFMM2D.
% Adjust the o.* options below to choose the experiment.
% Barnett 3/24/11
clear; addpath /home/alex/physics/leslie/gimbutas/lp2d/

N=1e3; k=10; verb=1;   % # qpadr pts, wavenumber, and verbosity.
s = segment.smoothstar(N,0.3,7); % actually has tight inside curvature, tricky!
%aj = zeros(1,50); bj = aj; aj(4) = -0.2; bj(end)=0.2;
%s = segment.smoothfourier(N,aj,bj);
de = domain([], [], s, -1);   % exterior
de.addlayerpot(s, [-1i*k 1], struct('quad','a','ord',8)); % CFIE, Alpert order
s.setbc(1, 'd', []);
p = scattering(de, []); p.setoverallwavenumber(k);
p.setincidentwave(pi/2 - pi/20);  % if just angle given, it's a plane wave

fprintf('testing N=%d; please wait about %g min...\n', N, N/1e4); 
o.FMM = 0; o.meth = 'iter'; %o.meth = 'direct';
tic; p.solvecoeffs(o); fprintf('solve done in %.3g sec\n', toc)
p.pointsolution(pointset(-1.5-1.5i))
if verb,  figure; tic; p.showthreefields(o); % Note FMM also used here!
  fprintf('solution field eval in %.3g sec\n', toc), end

