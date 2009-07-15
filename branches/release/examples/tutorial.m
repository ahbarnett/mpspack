% Example codes from tutorial document, also generates EPS figures for this doc
% Barnett 7/14/09

cd ../doc                  % so figures write out write tutorial.tex is
clear all classes
verb = 0;                  % if verb>0, generates EPS figures

% ------------- 2.. Laplace BVP ----------------

s = segment([], [0 1 0 2*pi]);
d = domain(s, +1);
d.k = 0;
d.addregfbbasis([], 8);
%f = @(z) real(exp(z)); % gives 1/n! coeffs in the real part
f = @(z) log(abs(z-2-3i));              % boundary data function of z=(x,y)
%d.k = 1; f = @(z) bessely(0,abs(z-1-2i));    % Helmholtz
s.setbc(-1, 'd', [], @(t) f(s.Z(t)));
% f = @(x,y) exp(x).*cos(y); s.setbc(-1, 'd', [], f); % another way to pass in f
p = bvp(d);
p.solvecoeffs;
p.bcresidualnorm
p.showsolution;
figure; opts.comparefunc = f; p.showsolution(opts);

if verb, % generate f:sd
  figure; set(gca,'fontsize', 14); s.plot;
  print -depsc2 seg.eps
  figure; set(gca,'fontsize', 14); opts.gridinside=0.05; d.plot(opts); axis off
  print -depsc2 dom.eps
  % generate f:u
  figure; set(gca,'fontsize', 14); p.showsolution; axis off;
  h=colorbar; set(h,'fontsize',20);
  print -depsc2 -painters u.eps
  figure; set(gca,'fontsize', 14); p.showsolution(opts); axis off
  h=colorbar; set(h,'fontsize',20);
  print -depsc2 -painters uerr.eps
  % problem: -painters stops the transparency being correctly rendered.
  close all
end


% --------------- 3.. Accuracy and convergence, radial domains ----------------
s.requadrature(50); p.solvecoeffs; p.bcresidualnorm

for N=1:15,
  d.bas{1}.N = N; p.solvecoeffs; r(N) = p.bcresidualnorm;
end
figure; semilogy(r, '+-'); xlabel('N'); ylabel('bdry err norm');

if verb, % generate f:conv
  figure; set(gca,'fontsize', 20); semilogy(r, '+-');
  xlabel('N'); ylabel('bdry err norm'); print -depsc2 N.eps
end

% radial function star-shaped domain
s = segment.radialfunc(50, {@(q) 1 + 0.3*cos(3*q), @(q) -3*0.3*sin(3*q)});
d = domain(s, +1);
d.k = 0;
d.addregfbbasis([], 8);
s.setbc(-1, 'd', [], @(t) f(s.Z(t)));
% f = @(x,y) exp(x).*cos(y); s.setbc(-1, 'd', [], f); % another way to pass in f
p = bvp(d);
p.solvecoeffs; p.bcresidualnorm
figure; opts.comparefunc = f; p.showsolution(opts);
if verb % generate f:radfunc
  h=colorbar; set(h,'fontsize',20); hold on; s.plot; axis off;
  print -depsc2 -painters radfunc.eps
end

% general analytic domain: crescent
a = 0.2; b = 0.8; w = @(t) exp(2i*pi*t);
s = segment(100, {@(t) w(t)-a./(w(t)+b), ...
                  @(t) 2i*pi*w(t).*(1 + a./(w(t) + b).^2)}, 'p');
figure; s.plot;  % kill the arrow maybe



% ------------- 4. Helmholtz, exterior -----------------------------------





% end
