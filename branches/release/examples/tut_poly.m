% Example codes from MPSpack tutorial, also generates EPS figures for this doc
% SECTION 5: POLYGONS

cd ../doc                             % so figures write out to doc/ directory
clear all classes; verb = 0;          % if verb>0, generates EPS figures

% triangle interior...
s = segment.polyseglist([], [1, 1i, exp(4i*pi/3)]);
tri = domain(s, 1);
if verb  % generate f:doms b
  figure; opts.gridinside=0.05; tri.plot(opts); axis off; print -depsc2 tri.eps
end

s.disconnect;                      % play with some domains
exttri = domain([], [], s, -1]);
s.disconnect; 
ss = s.translate(2);
exttwotri = domain([], [], {s(end:-1:1), ss(end:-1:1)}, {-1, -1});
%if verb  % generate f:doms c
%  figure; exttwotri.plot(opts); axis off; print -depsc2 exttwotri.eps
%end

% Helmholtz triangle BVP w/regfb...
clear all classes; verb = 0;                            % new START POINT
s = segment.polyseglist(50, [1, 1i, exp(4i*pi/3)]);
tri = domain(s, 1);
s.setbc(-1, 'd', [], @(t) 1+0*t);
tri.addregfbbasis(0, []); tri.bas{1}.rescale_rad = 1.0;
p = bvp(tri);
tri.k = 10;
Ns = 2:2:40; for i=1:numel(Ns)
  tri.bas{1}.N = Ns(i); p.solvecoeffs; r(i) = p.bcresidualnorm;
  %p.pointsolution(pointset(0)), nm(i) = norm(p.co); % check convergence
end
figure; loglog(2*Ns, r, '+-'); xlabel('# degrees of freedom'); ylabel('bdry err norm');
% figure; loglog(2*Ns, [r; nm], '+-');   % also plot coeff norm

if verb  % generate f:triconv a
  g = gcf;
  f = figure; fb = loglog(2*Ns, r, '+-'); set(gca,'fontsize', 20); axis tight;
  xlabel('# degrees of freedom'); ylabel('bdry err norm');
  print -depsc2 triFB.eps
  figure(g);
end

% Helmholtz triangle BVP w/nufb...          trying various basis corner combos
opts = []; opts.rescale_rad = 1.9; % err < 1e-8 is way too dependent on this!
opts.cornermultipliers = [1 0 0]; opts.type = 's';  % 'cs' not needed w/ reg FB
tri.addcornerbases([], opts);
Ns = 1:30; for i=1:numel(Ns)
  p.updateN(Ns(i)); nn(i) = p.N; % save the # dofs the problem used
  p.solvecoeffs; r(i) = p.bcresidualnorm; nm(i) = norm(p.co);
  p.pointsolution(pointset(0)) %converges to -3.9871841923
end
% add to previous figure...        (check norm w/ loglog(Ns, [r;nm], 'r+-');)
hold on; loglog(nn, r, 'r+-'); xlabel('N'); ylabel('bdry err norm');
% notice with two corner expansions there are 4N dofs in total

if verb, % generate f:triconv a and b
  figure(f); hold on; nufb = loglog(nn, r, 'r+-'); axis tight;
  legend([fb nufb], {'regular FB at origin', 'reg FB + corner nu-FB'},...
         'location', 'southwest');
  print -depsc2 triFB.eps
  figure; p.showsolution; axis off; h=colorbar; set(h,'fontsize',20);
  print -depsc2 -painters triu.eps
end

figure; tri.plot; tri.showbasesgeom;   % check corner bases correct

% above we notive value at z=0 converges to roughly same # correct digits as
% the boundary error norm.
