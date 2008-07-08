% test program for all basis types
% barnett 7/7/08
% adapted from ~/physics/leslie/latticesum/test_basis_and_multipole.m
%
% To Do:
% * prefer [A Ax Ay] output args
% * reg FB: fix |x|=0 NaN problem
% * check y-derivs too

dx = 0.01; g = -1:dx:1;              % plotting region
[xx yy] = meshgrid(g, g);
M = numel(xx);
p = pointset([xx(:) + 1i*yy(:)], ones(size(xx(:)))); % set up n for x-derivs

k = 10;

for type = 1:2;
  switch type
   case {1,2}             % ................. Reg FB
    N = 10;
    realflag = (type==1);
    b = regfbbasis(0, N, realflag, k);
    js = 1:b.Nf;             % indices to plot, here all of them
    fprintf('evaluating Reg FB basis... realflag=%d\n', realflag)
  end
  tic; [A An] = b.eval(p); t=toc;
  fprintf('\t%d evals in %.2g secs = %.2g us per eval\n',...
          numel(A), t, 1e6*t/numel(A))
  nnans = numel(find(isnan([u un])));
  if nnans, fprintf('\tproblem: # NaNs = %d\n', nnans); end
  n = numel(js);
  u = reshape(A(:,js), [size(xx) n]);   % make stack of rect arrays
  un = reshape(An(:,js), [size(xx) n]);        % from j'th cols of A, An
  nh = floor(sqrt(1.8*n)); nv = ceil(n/nh); % how many across and down, subplot
  
  % values...
  figure('name', sprintf('%d: Re[u_j(x)]', type)); c = .5;    % caxis range
  for j=1:numel(js), subplot(nv, nh, j); imagesc(g, g, real(squeeze(u(:,:,j))));
    caxis(c*[-1 1]); set(gca, 'ydir', 'normal'); axis equal tight, end
  subplotspace('vertical', -10); subplotspace('horizontal', -15);
  if ~realflag, figure('name', sprintf('%d: Im[u_j(x)]', type));
  for j=1:numel(js), subplot(nv, nh, j); imagesc(g, g, imag(squeeze(u(:,:,j))));
    caxis(c*[-1 1]); set(gca, 'ydir', 'normal'); axis equal tight, end
  subplotspace('vertical', -10); subplotspace('horizontal', -15);
  end
  
  % errors in x-derivs...   plots should have no values larger than colorscale
  unje = (u(:,3:end,:)-u(:,1:end-2,:))/(2*dx) - un(:,2:end-1,:); % FD error
  fprintf('\tmax err in u_x from FD calc = %g\n', max(abs(unje(:))))
  figure('name', sprintf('%d: FD err in Re[du_j/dn(x)]', type));
  c = c*(k*dx)^2; % typ err size of FD is (k.dx)^2 since 2nd order
  for j=1:numel(js), subplot(nv, nh, j);
    imagesc(g(2:end-1), g, real(squeeze(unje(:,:,j))));
    caxis(c*[-1 1]); set(gca, 'ydir', 'normal'); axis equal tight, end
  subplotspace('vertical', -10); subplotspace('horizontal', -15);
  if ~realflag, figure('name', sprintf('%d: FD err in Im[du_j/dn(x)]', type));
  for j=1:numel(js), subplot(nv, nh, j);
    imagesc(g(2:end-1), g, imag(squeeze(unje(:,:,j))));
    caxis(c*[-1 1]); set(gca, 'ydir', 'normal'); axis equal tight, end
  subplotspace('vertical', -10); subplotspace('horizontal', -15);
  end  
end