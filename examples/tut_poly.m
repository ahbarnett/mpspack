% Example codes from MPSpack tutorial, also generates EPS figures for this doc
% SECTION 5: POLYGONS

clear all classes; verb = 0;          % if verb>0, generates EPS figures

if 0 % triangle interior...
s = segment.polyseglist([], [1, exp(3i*pi/8), exp(5i*pi/4)]);
tri = domain(s, 1);
if verb  % generate f:doms b
  figure; opts.gridinside=0.05; tri.plot(opts); axis off;
  print -depsc2 ../doc/figs/tri.eps
end

s.disconnect;                      % play with some domains
exttri = domain([], [], s(end:-1:1), -1);
s.disconnect; 
ss = s.translate(2);
exttwotri = domain([], [], {s(end:-1:1), ss(end:-1:1)}, {-1, -1});
%if verb  % generate f:doms c
%  figure; exttwotri.plot(opts); axis off; print -depsc2 ../doc/figs/exttwotri.eps
%end
end

% Helmholtz triangle BVP w/regfb...
clear all classes; verb = 1;  % ============ possible code START POINT ========
%s = segment.polyseglist(100, [1, 1i, exp(4i*pi/3)]);
s = segment.polyseglist(50, [1, exp(3i*pi/8), exp(5i*pi/4)]);
tri = domain(s, 1);
s.setbc(-1, 'd', [], @(t) 1+0*t);
tri.addregfbbasis(0, []); tri.bas{1}.rescale_rad = 1.0;
p = bvp(tri); tri.k = 10;
Ns = 2:2:40; for i=1:numel(Ns)
  p.updateN(Ns(i)); p.solvecoeffs; r(i) = p.bcresidualnorm;
  p.pointsolution(pointset(0)), %nm(i) = norm(p.co); % check convergence
end
figure; loglog(2*Ns, r, '+-'); xlabel('# degrees of freedom'); ylabel('bdry err norm');
% figure; loglog(2*Ns, [r; nm], '+-');   % also plot coeff norm

if verb  % generate f:triconv a
  g = gcf;
  f = figure; fb = loglog(2*Ns, r, '+-'); set(gca,'fontsize', 20); axis tight;
  xlabel('# degrees of freedom'); ylabel('bdry err norm');
  print -depsc2 ../doc/figs/triFB.eps
  figure(g);
end

% Helmholtz triangle BVP w/nufb...          trying various basis corner combos
opts = []; opts.rescale_rad = 2.0; % err < 1e-8 is way too dependent on this!
%opts.cornermultipliers = [1 1 1];
tri.clearbases; tri.addcornerbases([], opts);
r = []; Ns = 1:13; for i=1:numel(Ns)
  p.updateN(Ns(i)); nn(i) = p.N; % save the # dofs the problem used
  p.solvecoeffs; r(i) = p.bcresidualnorm; %nm(i) = norm(p.co); % save coeff norm
  p.pointsolution(pointset(0))             % watch solution u(0)
end
% add to previous figure...        (check norm w/ loglog(nn, [r;nm], 'r+-');)
hold on; loglog(nn, r, 'r+-'); xlabel('N'); ylabel('bdry err norm');
% notice with two corner expansions there are 4N dofs in total

if verb, % generate f:triconv a and b
  figure(f); hold on; nufb = loglog(nn, r, 'r+-'); axis([4 100 1e-12 10]);
  legend([fb nufb], {'regular FB at origin', '3 corner nu-FB'},...
         'location', 'southwest');
  print -depsc2 ../doc/figs/triFB.eps
  figure; p.showsolution; axis off; h=colorbar; set(h,'fontsize',20);
  print -depsc2 -painters ../doc/figs/triu.eps
end

%figure; tri.plot; tri.showbasesgeom;   % check corner bases correct

% above we notive value at z=0 converges to roughly same # correct digits as
% the boundary error norm. However, we cannot get more than 11 digits on
% this value, -3.9871841923.
