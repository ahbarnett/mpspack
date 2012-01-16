% develop corner quadratures in MPSpack

if 0 % Prelim testing of the reparam funcs. Map (6.3) from Kress '91:
  s = 0:0.001:1; w=@(s) exp(-1./s)./(exp(-1./s)+exp(-1./(1-s)));
  figure; plot(s,w(s),'-');
  
  % map from Kress '91 which is more successful:
  q = 8; v=@(s) (1/q-1/2)*(1-s).^3 + (s-1)/q + 1/2; 
  w = @(s) v(s).^q ./ (v(s).^q + v(1-s).^q);
  s = 0:0.01:1; figure; plot(s,w(s), 'r-');
  figure; plot(s,v(s),'-'); hold on; plot(s,v(s).^q,'g-');
end

clear all classes;   % MPSpack segment 'pc' corner quad
N = 50; s = segment(N, [0 1], 'pc');
%figure; s.plot;

% Helmholtz BVP convergence... (interior or exterior, with pt src or const data)
o.kressq=4; s=segment.polyseglist(50, [1, exp(3i*pi/8), exp(5i*pi/4)], 'pc', o);
k = 100;      % choose wave#
inout = -1;  % choose expt:  +1 for exterior (nasty corners), -1 for interior
ftype = 0;   % choose bdry data type: ftype=0: f=H0; ftype=1: f=1
z0in = 0; z0out=1+1i; if inout==1, t=z0in; z0in=z0out; z0out=t; end % test pts
f = @(z) besselh(0,k*abs(z-z0out)); if ftype==1, f = @(z) 0*z + 1; end
if inout==-1, tri = domain(s, 1); else tri = domain([],[],s(end:-1:1),-1); end
for i=1:3, s(i).setbc(inout, 'd', [], @(t) f(s(i).Z(t))); end % f bdry data
tri.addlayerpot([], 'd');  % one layerpot per segment
p = bvp(tri); tri.k = k; Ns = 30:10:200;  % 30:10:150
u = nan*Ns;
for i=1:numel(Ns), N=Ns(i); p.updateN(N);      % convergence
  %min(abs(diff(vertcat(s.x)))) % too close?
  p.solvecoeffs; u(i) = p.pointsolution(pointset(z0in));
  fprintf('N=%d: u(0) = %.16g + %.16gi\n',N,real(u(i)),imag(u(i)))
end
if ftype==0, e = u-f(z0in); else   % decide ptwise error measure for u(z0in)
  if inout==-1, e = imag(u); else, e = u-u(end); end, end
figure; loglog(Ns,abs(e),'+-');hold on;plot(Ns,(Ns/10).^-4,'r-');axis tight;
plot(Ns,Ns.^-3,'m-');

%figure; imagesc(log10(abs(p.A))); colorbar  % show A
if 0, figure; o = []; if ftype==0, o.comparefunc = f; end, p.showsolution(o);
  p.plot; hold on; plot([z0in z0out], '+'); caxis(1e-12*[-1 1]);
end

if 0 % plot solution density two ways: wrt arclength, and wrt node index...
  ww = vertcat(s.w)'; figure; subplot(2,1,1);
  plot(cumsum(ww(:)), abs(p.co), '+-'); subplot(2,1,2); plot(abs(p.co),'+-');
end
