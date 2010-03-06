% test layerpot.evalfty by computing coeffs for y-slice lp rep in half-space,
% for some layerpot source on the other half-space.
% ............ relies on ftylayerpot.eval being correct
% barnett 2/24/10

clear all classes
src = 'd';      % 's' or 'd' for the source density type
rep = 's';      % 's' or 'd' for testing the half-space rep type by y-ray LP
side = -1;      % +-1, which side of x=0 the source is to be
om = 10;                                       % overall frequency
b = ftylayerpot(0, rep, struct('omega', om));  % make basis set to be tested
e = domain(); e.k = om; b.doms = e; % hook to an R^2 domain
s = segment([], [side*0.8 side*1.0+0.3i]);  % make horiz or vert to debug nx,ny
l = layerpot(s, src); l.doms = e; % to R^2
co = (6+24*(src=='s'))*ones(size(s.x));           % scale so u field is O(1)
dx = 0.03; x = -1:dx:1.5; y = -3:dx:3; [xx yy] = meshgrid(x,y);
p = pointset(xx(:)+1i*yy(:));
A = l.eval(p); u = reshape(A*co, size(xx)); % show only refl src
o=[]; %o.side = -sign(real(p.x));           % force the FTyLP evaluation side
[F Fx] = l.evalfty(b); E = b.eval(p,o);  % get the FTy u & u_n eval matrices
if rep=='d'
  ulp = -2*side*reshape(E*(F*co), size(xx)); % use FTy-DLP with tau = 2u
else ulp = 2*side*reshape(E*(Fx*co), size(xx)); end % or Fty-DLP, sig = -2u_n

figure; subplot(1,3,1); imagesc(x, y, real(u)); set(gca, 'ydir', 'normal');
axis equal tight; caxis([-1 1]); colorbar; b.showgeom;
xlabel('x'); ylabel('y');title(['LP src val: ' src]); colormap(jet(256));
subplot(1,3,2); imagesc(x, y, real(ulp)); set(gca, 'ydir', 'normal');
caxis([-1 1]); axis equal tight; colorbar; b.showgeom;
xlabel('x'); ylabel('y');title(['FTy val: ' rep]); colormap(jet(256));
subplot(1,3,3); imagesc(x, y, log10(abs(u-ulp))); set(gca, 'ydir', 'normal');
axis equal tight; xlabel('x');ylabel('y');title('log_{10} abs err');
caxis([-16 0]); colormap(jet(256)); colorbar;

%set(gcf, 'paperposition', [0 0 8 5]); print -depsc2 test_FTyDLP_refl.eps
% src=s,rep=d, sent to Leslie 2/24/10


