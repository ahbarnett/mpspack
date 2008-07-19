% BVP class test routine
% barnett 7/18/08, included exterior domain and Neumann 7/19/08

clear all classes
verb = 1; bctype = 'n';               % choose BC type: 'd' or 'n'
p = [1 1i exp(4i*pi/3)];              % triangle from Alex's inclusion paper
M = 40;
s(1) = segment(M, [p(1) p(2)]);
s(2) = segment(M, [p(2) p(3)], 't');
s(3) = segment(M, [p(1) p(3)]);       % note 3rd seg is backwards
pm = [1 1 -1]; d = domain(s, pm);     % build interior polygonal domain
k = 10; N = 30;
opts.real = 1; d.addrpwbasis(N, k, opts);    % choose a basis set for domain

xsing = 1 + 1i;                       % location of singularity
f = @(x) bessely(0, k*abs(x - xsing));  % an exact Helmholtz solution in domain
if bctype=='d'
  for j=1:numel(s)                      % Dirichlet: same func on all segments
    s(j).setbc(pm(j), 'd', [], f);      % note pm is needed: which side BC is on
  end
elseif bctype=='n'
  R = @(x) abs(x - xsing);
  fx = @(x) -k*real(x-xsing)./R(x).*bessely(1, k*R(x)); % components of grad f
  fy = @(x) -k*imag(x-xsing)./R(x).*bessely(1, k*R(x));
  for j=1:numel(s)                      % Neumann: same func on all segments
    gd = fx(s(j).x) .* real(s(j).nx) + fy(s(j).x) .* imag(s(j).nx);
    s(j).setbc(pm(j), 'n', [], gd);      % now pass data in rather than a func
  end

end
  
pr = bvp(d);
pr.solvecoeffs;
fprintf('L2 bdry error norm = %g, coeff norm = %g\n', pr.bcresidualnorm, ...
        norm(pr.co))
figure; subplot(1,2,1);
pr.showbdry; axis equal; title('BVP bdry and soln'); hold on; plot(xsing,'+');
opts.dx = 0.03; tic; [uN gx gy di] = pr.gridsolution(opts); toc
imagesc(gx, gy, uN); colorbar;set(gca,'ydir','normal'); subplot(1,2,2);
[xx yy] = meshgrid(gx, gy); zz = xx + 1i*yy;
fd = f(zz); ii = find(~isnan(di)); % ii is indices inside any domain
imagesc(gx, gy, uN - fd); colorbar; axis equal; set(gca,'ydir','normal');
title('soln err');
r = opts.dx * norm(uN(ii) - fd(ii));    % L2 interior error norm (est on grid)
fprintf('L2 interior error norm = %g\n', r);


if 0 % exterior domain ........................ 'd' case only
segment.disconnect(s); s = s(end:-1:1);      % flip segment order for exterior
pm = [1 -1 -1]; d = domain([], [], s, pm);   % new exterior polygonal domain
N = 10; opts.real = 1; d.addrpwbasis(N, k, opts);
f = @(x) besselj(0, k*abs(x));               % entire Helmholtz solution
for j=1:numel(s)                      % same func on all segments
  s(j).setbc(pm(j), 'd', [], f);      % note pm is needed: which side BC is on
end
pr = bvp(d);
pr.solvecoeffs;
fprintf('L2 bdry error norm = %g, coeff norm = %g\n', pr.bcresidualnorm, ...
        norm(pr.co))
figure; pr.showbdry; axis equal; title('BVP bdry'); hold on;
tic; [u gx gy di] = pr.gridsolution; toc, imagesc(gx, gy, u); caxis([-.5 .5]);
end

