% test program for all basis types: reg FB (incl fast), ...
% barnett 7/8/08
% adapted from ~/bdry/inclus/test_evalbasis.m.
%
% To Do:
% * prefer [A Ax Ay] output args ?
% * reg FB: fix |x|=0 NaN problem
% * check y-derivs too

clear all classes
verb = 1;
dx = 0.01; g = -1:dx:1;              % plotting region
[xx yy] = meshgrid(g, g);
M = numel(xx);
p = pointset([xx(:) + 1i*yy(:)], ones(size(xx(:)))); % set up n for x-derivs
k = 10;

for type = 1:3
  switch type
   case {1,2,3}             % ................. Reg FB: slow real/cmplx, fast
    N = 10;
    opts.real = (type==1); opts.fast = (type==3);
    b = regfbbasis(0, N, k, opts);
    js = 1:b.Nf;             % indices to plot, here all of them
    fprintf('evaluating Reg FB basis... real=%d, fast=%d\n', opts.real, opts.fast)
  end
  tic; [A An] = b.eval(p); t=toc;
  fprintf('\t%d evals in %.2g secs = %.2g us per eval\n',...
          numel(A), t, 1e6*t/numel(A))
  n = numel(js);
  u = reshape(A(:,js), [size(xx) n]);   % make stack of rect arrays
  un = reshape(An(:,js), [size(xx) n]);        % from j'th cols of A, An
  nnans = numel(find(isnan([u un])));
  if nnans, fprintf('\tproblem: # NaNs = %d\n', nnans); end
  nh = floor(sqrt(1.8*n)); nv = ceil(n/nh); % how many across and down, subplot
  
  if verb % values...
  figure('name', sprintf('%d: Re[u_j(x)]', type)); c = .5;    % caxis range
  for j=1:numel(js), subplot(nv, nh, j); imagesc(g, g, real(squeeze(u(:,:,j))));
    caxis(c*[-1 1]); set(gca, 'ydir', 'normal'); axis equal tight, end
  subplotspace('vertical', -10); subplotspace('horizontal', -15);
  if ~opts.real, figure('name', sprintf('%d: Im[u_j(x)]', type));
  for j=1:numel(js), subplot(nv, nh, j); imagesc(g, g, imag(squeeze(u(:,:,j))));
    caxis(c*[-1 1]); set(gca, 'ydir', 'normal'); axis equal tight, end
  subplotspace('vertical', -10); subplotspace('horizontal', -15);
  end
  end
  % errors in x-derivs...   plots should have no values larger than colorscale
  unje = (u(:,3:end,:)-u(:,1:end-2,:))/(2*dx) - un(:,2:end-1,:); % FD error
  fprintf('\tmax err in u_x from FD calc = %g\n', max(abs(unje(:))))
  if verb, figure('name', sprintf('%d: FD err in Re[du_j/dn(x)]', type));
  c = c*(k*dx)^2; % typ err size of FD is (k.dx)^2 since 2nd order
  for j=1:numel(js), subplot(nv, nh, j);
    imagesc(g(2:end-1), g, real(squeeze(unje(:,:,j))));
    caxis(c*[-1 1]); set(gca, 'ydir', 'normal'); axis equal tight, end
  subplotspace('vertical', -10); subplotspace('horizontal', -15);
  if ~opts.real, figure('name', sprintf('%d: FD err in Im[du_j/dn(x)]', type));
    for j=1:numel(js), subplot(nv, nh, j);
      imagesc(g(2:end-1), g, imag(squeeze(unje(:,:,j))));
      caxis(c*[-1 1]); set(gca, 'ydir', 'normal'); axis equal tight, end
    subplotspace('vertical', -10); subplotspace('horizontal', -15);
  end  
  end
end

