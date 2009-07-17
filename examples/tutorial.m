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

s.disconnect;
exttri = domain([], [], s, -1]);
s.disconnect; 
ss = s.translate(2);
exttwotri = domain([], [], {s(end:-1:1), ss(end:-1:1)}, {-1, -1});
%if verb  % generate f:doms c
%  figure; exttwotri.plot(opts); axis off; print -depsc2 exttwotri.eps
%end

% Helmholtz triangle BVP...
clear all classes
s = segment.polyseglist([], [1, 1i, exp(4i*pi/3)]); s.requadrature(50);
tri = domain(s, 1);
s.setbc(-1, 'd', [], @(t) 1+0*t);
tri.addregfbbasis(0, []); tri.bas{1}.rescale_rad = 1.0;
p = bvp(tri);
tri.k = 5;
Ns = 2:2:30; for i=1:numel(Ns)
  tri.bas{1}.N = Ns(i); p.solvecoeffs; r(i) = p.bcresidualnorm;
  %nm(i) = norm(p.co);
end
figure; loglog(Ns, r, '+-'); xlabel('N'); ylabel('bdry err norm');
%figure; loglog(Ns, nm, '+-'); xlabel('N'); ylabel('coeff norm');

if verb  % generate f:triconv a
  figure;loglog(Ns, r, '+-'); set(gca,'fontsize', 20); axis tight;
  xlabel('N'); ylabel('bdry err norm'); print -depsc2 triFB.eps
end




% end
