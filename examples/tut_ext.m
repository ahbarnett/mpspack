% Example codes from MPSpack tutorial, also generates EPS figures for this doc
% SECTION 4: EXTERIOR MFS

clear all classes; verb = 0;          % if verb>0, generates EPS figures
% exterior domain
tref = segment.radialfunc(100, {@(q) 1 + 0.3*cos(3*q), @(q) -0.9*sin(3*q)});
d = domain([], [], tref, -1);
d.k = 10; f = @(z) besselh(0,d.k * abs(z-0.3-0.2i));  % a radiative Helm soln
tref.setbc(1, 'd', [], @(t) f(tref.Z(t)));            % exterior Dirichlet data
% MFS
%d.addmfsbasis(tref.scale(0.8));  % rescaled segment copy
opts.tau = 0.06; d.addmfsbasis(tref, [], opts);
p = bvp(d);
% convergence plot
for N=5:5:80,
  p.updateN(N); p.solvecoeffs; r(N) = p.bcresidualnorm;
end
figure; semilogy(r, '+-'); xlabel('N'); ylabel('bdry err norm');

opts.comparefunc = f; figure; p.showsolution(opts);

if verb
  figure; set(gca,'fontsize', 20);
  semilogy(r, '+-'); xlabel('N'); ylabel('bdry err norm');
  print -depsc2 ../doc/extconv.eps
  figure; set(gca,'fontsize', 20);
  d.plot; p.showbasesgeom; p.showsolution; axis off;
  print -depsc2 ../doc/extgeom.eps
end

% multiply-connected domains. 1 hole...
tref.disconnect;                         % clears any domains from segment
c = segment([], [0.5 0.4 0 2*pi]);       % new circular segment
d = domain(tref, 1, c, -1);
% 2 holes...
tref.disconnect; c.disconnect;
smtref = tref.scale(0.3);                % create new rescaled copy of tref
smtref.translate(-0.3+0.5i);             % move the segment smtref
d = domain(tref, 1, {c smtref}, {-1 -1});
if verb  % generate f:doms a
  figure; set(gca, 'fontsize', 14); d.plot; axis off;
  print -depsc2 ../doc/twoholes.eps
end
  
% maybe solve MFS in multiply-connected case? Timo...
