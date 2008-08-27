% test basis.evallocalcopies, qpunitcell.evalbasescopies
% barnett 8/22/08

clear classes
k = 10; bas = ''; o.nei = 1;   % test params: bas = 'm', 'qL', 'qB' ('' to skip)

uc = qpunitcell(1, 0.5+1i, k, 10);  uc.setbloch(-1,-1);  % problem is pre-stored data from qpuclp lags behind uc.a in the naive summation.
if ~isempty(bas)
if bas=='m'
  N=10; t=2*pi*(1:N)/N; b = mfsbasis(0.2*exp(1i*t), [], [], k); % mfs basis
  g = -1.8:0.03:1.8;
elseif bas(1)=='q'
  b = qpuclayerpot(uc, bas(2), 's', k);  % both qL & qB tested ok
  g = -2.5:0.05:2.5;
end
[xx yy] = meshgrid(g,g); p = pointset(xx(:)+1i*yy(:), ones(size(xx(:))));
%u = b.eval(p);
tic; [A o.data] = b.evallocalcopies(p, uc, o);
fprintf('eval copies from cold = %.3g s\n', toc)
tic; A = b.evallocalcopies(p, uc, o);
fprintf('eval copies using stored data = %.3g s\n', toc)
%profile clear; profile on; for i=1:100, A = b.evallocalcopies(p, uc, o); end; profile off; profile viewer % test that matrix mults dominate reusing storage.
u = reshape(A*ones(size(A,2),1), size(xx));
figure; imagesc(g, g, real(u)); set(gca, 'ydir', 'normal'); colorbar;
hold on; uc.plot;
end

% now test the all-uc-bases code for block row of B matrix values...
uc.clearbases; uc.addqpuclayerpots;
g = -2.5:0.05:2.5;
[xx yy] = meshgrid(g,g); p = pointset(xx(:)+1i*yy(:), ones(size(xx(:))));
tic; [Br o.data] = uc.evalbasescopies(p, o);
fprintf('uc.evalbasescopies from cold = %.3g s\n', toc)
tic; Br = uc.evalbasescopies(p, o);
fprintf('uc.evalbasescopies using stored data = %.3g s\n', toc)
%profile clear; profile on; for i=1:10, Br = uc.evalbasescopies(p, o); end; profile off; profile viewer % test that matrix mults dominate reusing storage.
u = reshape(Br*ones(size(Br,2),1), size(xx));
figure; imagesc(g, g, real(u)); set(gca, 'ydir', 'normal'); colorbar;
hold on; uc.plot;

% test alpha-poly produces something...
M = 100; p = pointset(rand(M,1)+1i*rand(M,1), ones(M,1));
o.poly = 1; tic; [Br o.data] = uc.evalbasescopies(p, o);
fprintf('uc.evalbasescopies from cold, alpha-poly = %.3g s\n', toc)
tic; Br = uc.evalbasescopies(p, o);
fprintf('uc.evalbasescopies using stored data, alpha-poly = %.3g s\n', toc)
figure; showmatpoly(Br);

if 0 % interesting timing about stacking in block cols... better to preallocate
  tic; B = []; for i=1:10, B = [B ones(300,1000)]; end; toc  % 3-4 times slower!
  tic; B = zeros(300,10000); for i=1:10, B(:,1000*(i-1)+(1:1000))= ones(300,1000); end; toc
end

if 0  % 0.8 ms for setupbasisdofs
  tic; for i=1:100, [N noff] = uc.setupbasisdofs; end; toc
end