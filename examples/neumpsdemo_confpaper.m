% Demo Neumann MPS for eigenvalues. Incomplete.
% Barnett 6/15/10, adapted Specgeom PSPM conference paper w/ Hassell Feb 2011,
% adapted into MPSpack tree, 3/19/12

clear;
a = 0.3; b = 0.2; w = 3; % shape params (default for BNDS: a=.3 b=.2 w=3)
s = segment.radialfunc(450, {@(q) 1 + a*cos(w*(q+b*cos(q))), ...
                    @(q) -a*sin(w*(q+b*cos(q))).*w.*(1-b*sin(q)), ...
                    @(q) -a*cos(w*(q+b*cos(q))).*w^2.*(1-b*sin(q)).^2 + a*sin(w*(q+b*cos(q))).*w.*b.*cos(q)});  % includes curvature. 200
%s = segment.smoothstar(100,0,3); % circle
M = numel(s.w); L = sum(s.w); % # pts, perimeter
arcl = real(ifft(-1i*fft(s.speed).*[0 1./(1:M/2-1) 0 -1./(M/2-1:-1:1)]'));
arcl = arcl'/2/pi + L*s.t'; % spectral approx to sampled cumulative arclength
%figure; plot(s.t, s.w, '.-'); hold on; plot(s.t, arcl, 'g+-'); plot(s.t, cumsum(s.w), 'r-'); figure; plot(s.t, arcl'-cumsum(s.w), '-'); figure; plot(diff([arcl; arcl],3)); % cf crude O(1/M) arcl, check spectrally smooth
d = domain(s,1);
d.addmfsbasis(s, [], struct('tau',-0.025)); d.bas{1}.realflag=1; % for speed
p = bvp(d); p.updateN(100);  % 300 for paper
N = d.bas{1}.Nf; % figure; d.plot; d.showbasesgeom;
if 0, E = 2000; E*d.area/(4*pi) - d.perim/(4*pi)*sqrt(E) % est # via Weyl.
  xn = real(conj(s.x) .* s.nx);
  sxx = repmat(s.x,[1 numel(s.x)]); max(max(abs(sxx - sxx.'))) % diam of Omega
  M=400; M/(sqrt(E)*d.perim/2/pi), end       % ppw on bdry

%to = struct('meth','gt','eps',1e-14);  % gsvd, timo, Dirichlet-approx int norm
% NEUMANN CASE: 1 for naive, 2 for F_k corrected:
to = struct('meth','at','gep','e','eps',1e-15, 'neu', 1, 'arcl', arcl); % inf nrm, timo-reg, GEVP
%to = struct('meth','at','gep','v','eps',1e-14, 'neu', 1, 'arcl', arcl); % inf nrm, timo-reg, GEVP
%to = struct('meth','at','gep','e','eps',1e-14); % inf nrm, timo-reg, GEVP

if 0,%Es = 2000:.4:2100;%Es=50:.2:170; Es = 2096.240169+1e-0*[-5:5]; % ... sweep
  Es = 5:0.5:50; %Es = 8.3899:0.00001:8.3900; % 2nd min @ 8.38996
  ts = nan(N,numel(Es)); dets = nan*Es;
  for i=1:numel(Es), fprintf('E=%.16g\n', Es(i))
    t = sqrt(evp.tensionsq(d, Es(i), to)); ts(1:numel(t),i)=t; end

  figure; plot(Es, ts, 'k.', 'markersize', 5); hold on;
  plot(Es, min(ts,[],1), 'b.-');
  hold on; for i=0:10, plot(Es, 10*abs(besselj(i-1,sqrt(Es))-(i./sqrt(Es)).*besselj(i,sqrt(Es))),'r-'); end % for circle case only.
                                                                                  %axis([min(Es) max(Es) 0 max(min(abs(ts),[],1))]);
  axis([min(Es) max(Es) 0 sqrt(max(Es))]);
  xlabel E; ylabel('$t[\tilde{u}_{min}]$', 'interpreter', 'latex');
  set(gcf,'paperposition', [0 0 8 2]);
end  %  ...

%v=axis; v(4)=5; axis(v); title('disc neu MPS w/ F_k^{-1}'); print -depsc2 disc_neu_MPS_Fk.eps
%v=axis; v(4)=8; axis(v); title('disc neu MPS naive'); print -depsc2 disc_neu_MPS_naive.eps
%save rfn_neu_naive.mat
if 0, v=axis;v(4)=8; axis(v);
figure; hist(abs(diff(min(ts,[],1))), 40) % test how equal the slopes are
end

% plot a single Neumann eigenfunction
%E = 8.38996;
E = 2096.240170; p.updateN(400);N = d.bas{1}.Nf; % t < 1e-6 at 
%i=5; fzero(@(E) besselj(i-1,sqrt(E))-i*bessel(i,sqrt(E))./sqrt(E), 41) % circ E
%41.1601334801531
%E = 41.1601334801531; p.updateN(100);N = d.bas{1}.Nf; % t < 1e-6 at N=400
to = struct('meth','at','gep','e','eps',1e-15, 'neu', 1, 'arcl', arcl);
% note 'at' has eigvec bug!
% note gep=v method fails here. Not sure why.
[t V] = evp.tensionsq(d, E, to); sqrt(min(abs(t)))
j = find(abs(t)==min(abs(t)))
%j = 100;
p.co = V(:,j); % note many t's negative!
% check it has zero neumann BCs:
[A An] = d.bas{1}.eval(pointset(s.x,s.nx)); u = A*p.co; un = An*p.co;
figure; plot(arcl,un,'+-');
figure; [u gx gy] = p.showsolution(struct('dx', 0.01)); %no FMM meth w/ MFS yet!
%figure; contourf(gx, gy, abs(u).^2);axis equal; axis tight;
figure; imagesc(gx, gy, abs(u).^2); set(gca,'ydir','normal');
colormap(1-gray(256)); caxis([0 7e-5])
xlabel('x_1'); ylabel('x_2'); axis equal; axis tight;
hold on; x = [s.x; s.x(1)]; plot(real(x), imag(x), 'b-');
%title('E = \mu^2 = 2096.24016, tension t[u] = 3\times 10^{-6}');
set(gcf,'paperposition', [0 0 7 5]);
%print -depsc2 rfn_E2096.24016_mode_bw.eps


if 0,
  for E=2096.240165:1e-6:2096.240175, % sweep: it jumps around once <1e-6
    t = evp.tensionsq(d, E, to); j = find(abs(t)==min(abs(t))); t(j)
  end
end

