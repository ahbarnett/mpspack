% test Fourier-transform y Sommerfeld-integral layer pot
% barnett 2/23/10

clear all classes
om = 10;                                       % overall frequency
b = ftylayerpot(0, 's', struct('omega', om));  % make basis set to be tested
e = domain(); e.k = om; b.doms = e; % hook to an R^2 domain
%b.updateN(150);   % bump up according to various omega (eg <<1 or >>1)

if 0 % evaluation error test, using analytic derivatives...
m = mfsbasis(pointset(0, 1)); m.doms = e; m.eta = Inf; % choose S/D (eta=Inf,0)
dx = 0.03; x = -2:dx:-.5; %.5:dx:2;                    % choose + or - for x
y = -2:dx:2; [xx yy] = meshgrid(x,y);                  % choose grid
p = pointset(xx(:)+1i*yy(:), exp(1i*2*pi*rand(size(xx(:))))); % random normals 
[Am Amx Amy] = m.eval(p); [Am Amn] = m.eval(p); % exact point src
if m.eta==0, b.a=[0 1]; end, [A Ax Ay] = b.eval(p); [A An] = b.eval(p);% eta mix
co = ones(numel(b.kj),1) / 2/pi;        % coeffs giving i/4 H_0 pt src in SLP
nam={'u', 'u_x', 'u_y', 'u_n'}; for i=1:4       % loop over u, x-, y-, n-deriv
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

if 1, t = qpstrip(1, om); %figure; t.plot; % test qpstrip & Q matrix...
  t.addqpftylayerpots;

% test Q matrix --------------------------------------------------------------
t.setbloch(exp(1i*pi/7));    % random Bloch, will dep on u_inc plane wavevector
Q = t.evalbasesdiscrep(); min(svd(Q)) % test Q goes singular where predicted...
t.setbloch(exp(1i*om*(real(t.e)))); Q = t.evalbasesdiscrep(); min(svd(Q))
if 0, ps = -pi:0.1:pi; ss = nan*ps;  % sweep Bloch phase & plot min sing val
  for i=1:numel(ps), t.setbloch(exp(1i*ps(i)));
    Q = t.evalbasesdiscrep(); ss(i) = min(svd(Q)); end
    figure; plot(ps, ss, '+-'); vline(mod(om*(real(t.e))+pi, 2*pi)-pi);
end
end

% set mpspack up so once qpstrip defined, can call a scattering bvp and solve
