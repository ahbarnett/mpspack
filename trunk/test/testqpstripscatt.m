% test toy Dirichlet obstacle scatt by hand (no bvp class yet), BWLP ikS+D on pO
% barnett 2/26/10. tried wobbly Dirichlet wall, 3xUC, 3/3/10

clear all classes; v = 0;        % verbosity
%thi = -pi/5; om = 10; %7.76644415490187; %(t.a=1) % inc ang, overall freq
%thi = -pi/3; om = 4*pi; % 'double' Woods anomaly, with Bloch a=1
om = 10; thi = -acos(1 - 2*pi/om); % 'single' Wood's anomaly (generic Bloch a)
N = 80; M = 300;                 % N on obst LP, M on FTy LPs (M>300 for Wood)
o.nei = 2; o.buf = 0;            % scheme params: direct sum dist, strip buffer
kvec = om*exp(1i*thi);
ui = @(x) exp(1i*real(conj(kvec) .* x)); %uix = @(x) 1i*real(kvec)*ui(x);
d=1; t = qpstrip(d*(1+2*o.buf), om); % d = x-periodicity (not t.e for buf>0)
a = exp(1i*real(conj(kvec) * d)); t.setbloch(exp(1i*real(conj(kvec) * t.e)));
s = scale(segment.smoothstar(N, 0.3, 3), 0.25); % closed obst (0.39978=Neu res)
%h=.2; s = segment(N, {@(t) .5-t+1i*h*sin(2*pi*t), @(t) -1+2i*pi*h*cos(2*pi*t), @(t) -4i*pi^2*h*sin(2*pi*t)},'p'); % wobbly wall (note R-to-L to get normal up)
% allow translation, and mult uinc by phase so hits obst at same phase...
tr = 0; %0.02; s.translate(tr); ui = @(x) ui(x)*exp(-1i*tr*real(kvec));
de = domain([], [], s, -1);   % infinite domain exterior to inclusion
%de.k = om; de.addlayerpot(s, [i*om 1]); l=de.bas{1}; %abbrev for obst layerpot
de.k = om; de.addlayerpot(s, 'd'); l=de.bas{1}; %abbrev for surf layerpot
t.addqpftylayerpots(struct('M',M)); % override default M

if v, dx = 0.025; x = -1.5+1e-3:dx:1.5; y = -2:dx:2;  % plot grid (shifted)
[xx yy] = meshgrid(x,y); p = pointset(xx(:)+1i*yy(:));
uig = reshape(ui(p.x),size(xx)); figure; title('u_{inc}');
set(gcf,'position',get(gcf,'position') + [0 0 300 0]); subplot(1,3,1);
imagesc(x, y, real(uig)); hold on; s.plot; t.plot; colormap(jet(256));
set(gca,'ydir','normal'); axis equal tight; caxis([-2 2]); hold on;
plot([0 -cos(thi)],[0 -sin(thi)], 'k-', 'linewidth', 3);  % show inc direc
n = 100; for i=-n:n, cn = cos(thi)+i*2*pi/d/om; % cos
  if abs(cn)<=1,plot([0 cn],[0 -sqrt(1-cn^2)],'m-','linewidth',3);end,end%Bragg
title('u_{inc} & Bragg directions'); end

tic; % set up matrix blocks... for direct copies translate targets always
A = l.eval(s, struct('dom', de));  % needs domain to eval in (selfint)
for n=[-o.nei:-1 1:o.nei], A = A + a^n * l.eval(pointset(s.x-d*n)); end
% note B is not in the [tau; -sigma] hat order, rather bas order...
B = [t.bas{1}.eval(s) + t.a*t.bas{1}.eval(pointset(s.x-t.e)), ...        % Shat
     t.bas{2}.eval(s) + t.a*t.bas{2}.eval(pointset(s.x-t.e))];          % Dhat
C = zeros(2*M,N); f = utils.copy(t.bas{1}); % dummy FTyLP for translated evals
for n=o.nei-2*o.buf:o.nei   % general o.buf case (o.nei>=o.buf)
  f.orig = t.Lo-n*d; [Cj Cxj] = l.evalfty(f); C = C + a^n*[Cj;Cxj];
  f.orig = t.Lo-(-1-2*o.buf-n)*d; [Cj Cxj] = l.evalfty(f);
  C = C - a^(-1-2*o.buf-n)*[Cj;Cxj]; end  % careful, since L's phase = 1 always
Q = t.evalbasesdiscrep();
E = [A B; C Q]; fprintf('E filled in %g s\n', toc); 
%fprintf('||E||_1=%g, cond(E)=%g, min sing val Q = %g\n', ...
%        norm(E,1), cond(E), min(svd(Q)));
tic; co = E \ [-ui(s.x); zeros(size(Q,1),1)];      % RHS is [-u_inc; 0]
fprintf('linear solve done in %g s\n', toc); 
%co = rand(N+2*M, 1);           % dummy for now
%figure; plot(abs([co(1:N); co(1:N)]), '+-'); % test periodicity of wobbly tau

%figure; plot([real(co) imag(co)], '+-')
%title('\omega=10.6 quad=u showing 1/sqrt sings at each Bragg k_y!')


if v, Ag = l.eval(p); % Plot soln: evaluation matrix for obst copies + FTyLPs
for n=[-o.nei:-1 1:o.nei], Ag = Ag + a^n * l.eval(pointset(p.x-d*n)); end
Lg = [t.bas{1}.eval(p) + t.a*t.bas{1}.eval(pointset(p.x-t.e)), ...  %FTyLP eval
      t.bas{2}.eval(p) + t.a*t.bas{2}.eval(pointset(p.x-t.e))];   % buf=0
ug = reshape(Ag*co(1:N) + Lg*co(N+1:end), size(xx));  % compute u_s
subplot(1,3,2); imagesc(x, y, real(ug)); hold on; s.plot; t.plot; title('u_s');
colormap(jet(256)); set(gca,'ydir','normal'); axis equal tight; caxis([-2 2]);
if o.buf==0, % wrap soln correctly for |x|<1.5... needs dx to divide width (e)
  ep=0; for i=[-1 1], ug(find(xx>=t.Lo+i*d+ep&xx<=t.Ro+i*d+ep)) = ...
    ug(find(xx>=t.Lo+ep&xx<=t.Ro+ep)) * a^i; end, end
subplot(1,3,3); imagesc(x, y, real(ug+uig)); title('u total wrapped');
hold on; s.translate(-d).plot; s.plot; s.translate(d).plot;
colormap(jet(256)); set(gca,'ydir','normal'); axis equal tight; caxis([-2 2]);
%print -depsc2 s.14_a1_w7.76_first.eps
end  % note field inside obst, or below the wobbly wall, is plotted wrong.

loc = 0.3+0.8i; p = pointset(loc+tr); Ag = l.eval(p); % val at single test point
for n=[-o.nei:-1 1:o.nei], Ag = Ag + a^n * l.eval(pointset(p.x-d*n)); end
Lg = [t.bas{1}.eval(p) + t.a*t.bas{1}.eval(pointset(p.x-t.e)), ...  %FTyLP eval
      t.bas{2}.eval(p) + t.a*t.bas{2}.eval(pointset(p.x-t.e))];   % buf=0
u = Ag*co(1:N) + Lg*co(N+1:end);  % compute u_s
fprintf('u_s (tau+FTyLP) at (%.3g,%.3g) = %.16g + %.16g i\n', ...
        real(loc),imag(loc), real(u), imag(u));

y=-1:.01:1; % Test: eval on L,R themselves: use opts.side in ftylp.eval ...
pt = pointset([t.Lo+1i*y'; t.Ro+1i*y'], 1+0*[y y]'); % L & R, rightwards normal
% it's getting tedious to keep summing this way...
[At Ant] = l.eval(pt); % evaluation w/ derivs, for obst copies + FTyLPs
for n=[-o.nei:-1 1:o.nei], [Aj Anj] =  l.eval(pointset(pt.x-d*n, pt.nx));
  At = At + a^n * Aj; Ant = Ant + a^n * Anj; end
o.side = [ones(size(y')); ones(size(y'))];  % force signs for L, L+e jump rels
[Lj1 Lnj1] = t.bas{1}.eval(pt,o); [Lj2 Lnj2] = t.bas{2}.eval(pt,o); % u_n too!
Lt = [Lj1 Lj2]; Lnt = [Lnj1 Lnj2];  % sources from L
o.side = [-ones(size(y')); -ones(size(y'))];  % jump rels, everything to left
[Lj1 Lnj1] = t.bas{1}.eval(pointset(pt.x-t.e, pt.nx),o); % srcs from R (need nx)
[Lj2 Lnj2] = t.bas{2}.eval(pointset(pt.x-t.e, pt.nx),o);
Lt = Lt + t.a*[Lj1 Lj2]; Lnt = Lnt + t.a*[Lnj1 Lnj2]; 
ut = At*co(1:N) + Lt*co(N+1:end); uL = ut(1:numel(y)); uR = ut(numel(y)+1:end);
unt = Ant*co(1:N) + Lnt*co(N+1:end);unL=unt(1:numel(y));unR=unt(numel(y)+1:end);
if v, figure; semilogy(y, [abs(uL - uR/t.a) abs(unL - unR/t.a)],'+-');
title('QP errors in u_s'); xlabel('y'); legend('discrep f','discrep f'''); end
fprintf('max QP u err = %g, u_n err = %g\n', max(abs(uL - uR/t.a)), ...
        max(abs(unL - unR/t.a)))
%figure; subplot(1,2,1); plot(y, real([uL uR]), '+-'); % test slice fields
%subplot(1,2,2); plot(y, imag([uL uR]), '+-');

%set(gcf,'paperposition', [0 0 8 5]); print -depsc2 qpstrip_scatt_w10.eps

if v>1   % check conditioning of E: r sing vecs with small sing vals are due to
         % Q, around k=0 on the Sommerfeld contour...
  figure; imagesc(real(E)); colormap(jet(256)); caxis(.1*[-1 1]); title('E');
  [U S V] = svd(E);
  figure; imagesc(abs(V)); title('V');
  figure; semilogy(abs(V(:,end))); title('smallest r sing vec magnitude');
end
