% Example codes from MPSpack tutorial, also generates EPS figures for this doc
% SECTION 3: CONVERGENCE

cd ../doc                             % so figures write out to doc/ directory
clear all classes; verb = 0;          % if verb>0, generates EPS figures

% code from previous section used for set-up...
s = segment([], [0 1 0 2*pi]);
d = domain(s, +1);
d.k = 0;
d.addregfbbasis([], 8);
f = @(z) log(abs(z-2-3i));              % boundary data function of z=(x,y)
s.setbc(-1, 'd', [], @(t) f(s.Z(t)));
p = bvp(d);

% new code
s.requadrature(50); p.solvecoeffs; p.bcresidualnorm

% convergence plot
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

