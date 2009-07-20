% Example codes from MPSpack tutorial, also generates EPS figures for this doc
% SECTION 2: LAPLACE

cd ../doc                             % so figures write out to doc/ directory
clear all classes; verb = 0;          % if verb>0, generates EPS figures

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

