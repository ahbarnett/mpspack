% test Fourier-transform y Sommerfeld-integral layer pot
% barnett 2/23/10, high-aspect added 4/28/10

clear all classes
om = 10; maxy = 5;                   % overall frequency, vertical extent
b = ftylayerpot(0, 's', struct('omega',om, 'maxy',maxy, 'minx',.7)); % estim's M
e = domain(); e.k = om; b.doms = e; % hook to an R^2 domain
%b.updateN(400);    % override M

if 0 % evaluation error test, using analytic derivatives...
m = mfsbasis(pointset(0, 1)); m.doms = e; m.eta = Inf; % choose S/D (eta=Inf,0)
dx = 0.03; x = -2:dx:-.5; %.5:dx:2;                    % choose + or - for x
y = -.5:dx:maxy; [xx yy] = meshgrid(x,y);                % choose grid, large y
p = pointset(xx(:)+1i*yy(:), exp(1i*2*pi*rand(size(xx(:))))); % random normals 
[Am Amx Amy] = m.eval(p); [Am Amn] = m.eval(p); % exact point src
if m.eta==0, b.a=[0 1]; end, [A Ax Ay] = b.eval(p); [A An] = b.eval(p);% eta mix
co = ones(numel(b.kj),1) / 2/pi;        % coeffs giving i/4 H_0 pt src in SLP
nam={'u', 'u_x', 'u_y', 'u_n'}; for i=1:2       % loop over u, x-, y-, n-deriv
  if i==1, u=A*co;uex=Am; elseif i==2, u=Ax*co;uex=Amx;  % get relevant field
  elseif i==3, u=Ay*co;uex=Amy; else u=An*co;uex=Amn; end
  uex = reshape(uex, size(xx)); u = reshape(u, size(xx));
  figure; set(gcf, 'name', sprintf('test eta=%g, %s', m.eta, nam{i}));
  subplot(1,2,1); imagesc(x, y, real(u)); set(gca, 'ydir', 'normal');
  xlabel('x');ylabel('y');title('FTySLP'); axis equal; colormap(jet(256));
  subplot(1,2,2); imagesc(x, y, log10(abs(u-uex))); set(gca, 'ydir', 'normal');
  xlabel('x'); ylabel('y'); title('log_{10} abs err'); caxis([-16 0]);
  colormap(jet(256)); colorbar; axis equal;
end
end

t = qpstrip(1, om); %figure; t.plot;
t.addqpftylayerpots;

if 0 % ... test qpstrip & Q matrix...
  disp('test the Q matrix (FtyLP discrepancy on a qpstrip):');
t.setbloch(exp(1i*pi/7));    % random Bloch, will dep on u_inc plane wavevector
Q = t.evalbasesdiscrep();
fprintf('min sing val Q, typical alpha (O(1)): %.3g\n', min(svd(Q)))
t.setbloch(exp(1i*om*(real(t.e)))); Q = t.evalbasesdiscrep();
fprintf('min sing val Q, Wood''s alpha (small): %.3g\n', min(svd(Q)))
if 0, ps = -pi:0.1:pi; ss = nan*ps;  % sweep Bloch phase & plot min sing val
  for i=1:numel(ps), t.setbloch(exp(1i*ps(i)));
    Q = t.evalbasesdiscrep(); ss(i) = min(svd(Q)); end
    figure; plot(ps, ss, '+-'); vline(mod(om*(real(t.e))+pi, 2*pi)-pi); end
end

if 1 % test qpftylayerpot by eval on a grid
dx = 0.03; x = -.5:dx:1; %.5:dx:2;                    % choose + or - for x
y = -.5:dx:maxy; [xx yy] = meshgrid(x,y);                % choose grid, large y
p = pointset(xx(:)+1i*yy(:));
A = t.bas{2}.eval(p); co = (1./(t.bas{1}.lp.kj-1i)).'; u = A*co;
u = reshape(u, size(xx));
figure; imagesc(x, y, real(u)); set(gca, 'ydir', 'normal'); caxis(10*[-1 1]);
xlabel('x');ylabel('y');title('QP FTySLP'); axis equal; colormap(jet(256));
t.bas{1}.showgeom;
end
