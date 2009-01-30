% test the spurious distance function used to try to kill spurious modes
% was for 3xUC scheme  - NOW OBSOLETE!
% barnett 1/9/09

clear all classes
k = 5;
dummyMu = 30;
uc = qpunitcell(1, 0.5+1i, k, dummyMu);       % decide the UC shape
%uc.spuriousdistance, stop

da = 0.01;
aangs = pi*(-1:da:1);    % alpha angles  (make 0 for single slice)
bangs = pi*(-1:da:1);     % beta angles
f = NaN*zeros(numel(aangs), numel(bangs));
for i=1:numel(aangs)
  for j=1:numel(bangs)
    uc.setbloch(exp(1i*aangs(i)), exp(1i*bangs(j)));
    f(i,j) = uc.spuriousdistance;
  end
end
figure; uc.setbloch(1, 1); uc.showbrillouin; hold on; R = uc.recip/(2*pi);
[aa bb] = meshgrid(aangs, bangs); ab = R*[aa(:)'; bb(:)'];
aa = reshape(ab(1,:), size(aa)); bb = reshape(ab(2,:), size(aa));
surf(aa, bb, 0*aa, f');
colorbar; shading interp; axis equal; set(gcf, 'name', 'Brillouin zone');
set(gca, 'ydir', 'normal'); xlabel('Bloch k_x'); ylabel('Bloch k_y');
v = axis;
for i=-2:1/3:2, for j=-3:1/3:3, x = i*uc.r1+j*uc.r2;
    if floor(i)-i~=0 | floor(j)-j~=0, circle(real(x),imag(x),k,'k-'); end
end, end              % spurious k-circles
axis(v);

