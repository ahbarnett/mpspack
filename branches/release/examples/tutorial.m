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

% exterior domain
tref = segment.radialfunc(50, {@(q) 1 + 0.3*cos(3*q), @(q) -0.9*sin(3*q)});
d = domain([], [], tref, -1);
   
% TIMO add code for mfs

%...

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
  print -depsc2 twoholes.eps
end
  
% ----------------------- 5. Corners -------------------------------------

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

% Helmholtz triangle BVP w/regfb...  (may start code afresh here)
clear all classes; verb = 0;
s = segment.polyseglist(50, [1, 1i, exp(4i*pi/3)]);
tri = domain(s, 1);
s.setbc(-1, 'd', [], @(t) 1+0*t);
tri.addregfbbasis(0, []); tri.bas{1}.rescale_rad = 1.0;
p = bvp(tri);
tri.k = 10;
Ns = 2:2:40; for i=1:numel(Ns)
  tri.bas{1}.N = Ns(i); p.solvecoeffs; r(i) = p.bcresidualnorm;
  p.pointsolution(pointset(0)), %nm(i) = norm(p.co); % check convergence
end
figure; loglog(2*Ns, r, '+-'); xlabel('# degrees of freedom'); ylabel('bdry err norm');

if verb  % generate f:triconv a
  g = gcf;
  f = figure; fb = loglog(2*Ns, r, '+-'); set(gca,'fontsize', 20); axis tight;
  xlabel('# degrees of freedom'); ylabel('bdry err norm');
  print -depsc2 triFB.eps
  figure(g);
end

% Helmholtz triangle BVP w/nufb...          trying various basis corner combos
tri.clearbases;
opts = []; opts.rescale_rad = 1.9; % err < 1e-8 is way too dependent on this!
opts.cornerflags = [1 0 0]; opts.type = 's';  % 'cs' not needed if reg FB used
tri.addcornerbases([], opts);
opts = []; opts.rescale_rad = 1.0; tri.addregfbbasis(0, [], opts); % reg FB too
Ns = 1:30; for i=1:numel(Ns)
  for j=1:numel(tri.bas), tri.bas{j}.N=Ns(i); end; nn(i) = p.N; % # dofs
  p.solvecoeffs; r(i) = p.bcresidualnorm; nm(i) = norm(p.co);
  p.pointsolution(pointset(0)) %converges to -3.9871841923
end
% add to previous figure...        (check norm w/ loglog(Ns, [r;nm], 'r+-');)
hold on; loglog(nn, r, 'r+-'); xlabel('N'); ylabel('bdry err norm');
% notice with two corner expansions there are 4N dofs in total

if verb, % generate f:triconv a and b
  figure(f); hold on; nufb = loglog(4*Ns, r, 'r+-'); axis tight;
  legend([fb nufb], {'regular FB at origin', 'fractional FB at corners'},...
         'location', 'southwest');
  print -depsc2 triFB.eps
  figure; p.showsolution; axis off; h=colorbar; set(h,'fontsize',20);
  print -depsc2 -painters triu.eps
end

figure; tri.plot; tri.showbasesgeom;   % check corner bases correct

% above we notive value at z=0 converges to roughly same # correct digits as
% the boundary error norm.

% end
