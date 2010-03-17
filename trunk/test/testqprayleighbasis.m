% test the Rayleigh basis in a half-strip (qpstrip)
% Barnett 3/9/10

clear all classes; v = 1;        % verbosity
thi = -pi/5; om = 10; %7.76644415490187; %(t.a=1) % inc ang, overall freq
%om = 10; thi = -acos(1 - 2*pi/om); % 'single' Wood's anomaly (generic Bloch a)

kvec = om*exp(1i*thi); 
d = 1; yT = 1; T = segment([], 1i*yT+[d/2,-d/2], 'g'); up = -1;
tt = qpstrip(d, om, struct('seg',T, 'pm',-up)); % semi-bounded strip (-1 is up)
tt.setbloch(exp(1i*real(conj(kvec) * tt.e))); % Bloch
a = exp(1i*real(conj(kvec) * d)); tt.addqprayleighbasis(4); b = tt.bas{1};
%figure; tt.plot(struct('gridinside', .1)); tt.showbasesgeom;

if v, gx = -d/2:.01:d/2; gy = yT + (-2:0.04:2);  % plotting region
[xx yy] = meshgrid(gx, gy);
p = pointset([xx(:) + 1i*yy(:)], 1i*ones(size(xx(:)))); % set up n for y-derivs
[A An] = b.eval(p);
figure;
for j=1:b.Nf, subplot(3,3,j);  % show a bunch of them
co = zeros(size(A,2),1); co(j) = 1; u = reshape(A*co, size(xx)); % eval field
imagesc(gx, gy, real(u)); caxis([-1 1]); set(gca, 'ydir', 'normal');
colormap(jet(256)); axis tight equal; tt.plot; hold on;
n=100; for i=-n:n, cn = cos(thi)+i*2*pi/d/om; % cos used to get Bragg angs
  if abs(cn)<=1,plot([0 cn],[0 up*sqrt(1-cn^2)],'m-','linewidth',3);end,end
end
end
