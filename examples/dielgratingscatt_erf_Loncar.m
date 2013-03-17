% Scattering by plane wave from z-invariant (2D) periodic dielectric interface.
% Corners and multiple segments (in s) to represent the interface. TM only
% Wood's case only for upper layer region. Based upon fig_dielinterface 2/5/11
% Uses MPSpack V1.32 or later, and hacks the qp classes for now. Barnett 3/16/13
%
% To the user: please ignore domain warnings spat out by MPSpack.
% Adjustable parameters: om (freq), inc. wave angle, d (period), surface shape,
%  v (0:no plot, 1:plot), refr (index of lower layer),
%  N (incr to reduce error), M (unlikely; incr to reduce "rad cond" err).
%
% To do: * why is Kress not working with 'pc' in this QP setting? Alpert does
%         (Alpert easier to code anyway, since doesn't need Zpp func)
% * currently the central set of segs must lie in |x|<=0.5 otherwise large errs.

clear; v = 1; d = 1.0; % verbosity (0,1,...), spatial period
N = 100; o.kressq=5; % default pts per seg, and Kress corner bunching param
resc = 1;  % whether to do row+col rescaling of lin sys.
%  ========== Choose interface shape via 1 or more segments:
%h=1*.4; vt=0; s = segment(N, {@(t) 1i*vt + (.5-t)*d+1i*h*sin(2*pi*t), @(t) -d+2i*pi*h*cos(2*pi*t), @(t) -4i*pi^2*h*sin(2*pi*t)},'pc',o); % note R-to-L sense
%
%h=0; s = segment(N,[.5*d+1i*h,-.5*d+1i*h],'pc',o);  % flat
%
%t=0.45*pi;R=d/2/sin(t);h=R*(1-cos(t)); % arc w/ top corner (t->pi/2 bad)
%s=segment(N,[1i*R R -pi/2+t -pi/2-t],'pc',o); %(cont)
%
%h=0.6; s = [segment(N,[.5*d,-.3*d+1i*h],'pc',o) segment(N,[-.3*d+1i*h,-.5*d],'pc',o)];  % flat or zigzag w/ 2 segments, must be in R to L order
%
%h=.5; s=[segment(N, [.25*d .25*d 0 -pi], 'pc',o) segment(N, [-.25*d .25*d 0 pi], 'pc',o)]; % glued semicircles, C^1 
%
if 1, % use general horiz filling fraction funcs of vertical coord y: (Loncar)
  h = 2.0; % h: stopping half-height, |y|<h/2
  a=0.6; f = @(y) (1-erf(y/a))/2; % a vert scale, erf (deriv has Fourier decay)
  fp = @(y) (-1/sqrt(pi)/a)*exp(-(y/a).^2); % its deriv. Gauss approx to Dolph
  %f = @(y) 0.5-y/h; fp = @(y) -1/h+0*y; % tri wave to test
  %f = @(y) 0.5+0*y; fp = @(y) 0*y; % rounded vertical walls
  %a=0.5; b=1/erf(h/2/a); % expand by factor b to make a true non-cusp corner
  %f = @(y) (1-b*erf(y/a))/2; % a vert scale, f fill frac in [0,1] vs y
  %fp = @(y) b*(-1/sqrt(pi)/a)*exp(-(y/a).^2);
  N = 150; o.kressq=3;
  s(1) = segment(N, {@(t) .5*d-.5*d*f(h*(.5-t))+1i*h*(.5-t), @(t) h*.5*d*fp(h*(.5-t))-1i*h, @(t) 0*t}, 'pc',o); % dummy curvature (Alpert doesn't need)
  s(2) = segment(N, {@(t) -.5*d+.5*d*f(h*(t-.5))+1i*h*(t-.5), @(t) h*.5*d*fp(h*(t-.5))+1i*h, @(t) 0*t}, 'pc',o);
  r = f(h/2)*d/2; fprintf('f spike rounding radius = %.3g\n',r)
  if r>1e-12, s(3)=s(2); % if there's a gap in the function, round
    o.kressq=3; Nr=30; s(2) = segment(Nr, [-.5i*h r 0 -pi], 'pc',o);
    s(4) = segment(Nr, [-.5*d+.5i*h r 0 pi], 'pc',o); % rounding -> 4 segs tot
    s = s.translate(r); end % this keeps all x coords in (-.5,.5)
end
% =========== (done choosing shape)
if v>1, figure; s.plot; drawnow, end
fresnel = (h==0); % if flat interface, debug by showing error from analytic
om = 5; % freq (<20, say), overall wavenumber = 2pi/lambda in upper material
nei = 1; buf = 0; M = 110;  % don't expect buf to work yet! M: periodizing dofs

du = domain([], [], s(end:-1:1), -1); du.cloc=0; % (most of) R2 in this domain
dd = domain([], [], s, 1); dd.cloc=0; % (most of) R2 in this domain
refr = 2;         % refractive index of lower layer
dd.setrefractiveindex(refr); omd = om*refr;
s.setmatch('diel', 'TM');
o.quad = 'a'; o.ord = 16;       % Alpert near-diagonal correction
%o.quad = 'm';                 % Kress correction, doesn't work
s.addinoutlayerpots('d', o);                      % double-sided layerpot
s.addinoutlayerpots('s', o);                      % "

p = scattering(du, dd); dd.isair=2; % scatt prob used for A matrix; new isair!
p.setoverallwavenumber(om);
p.setincidentwave(-pi/3);   % incident wave angle from horizontal, in (-pi,0)
p.fillquadwei;             % this is only for obstacle mismatch (blocks A,B)
p.sqrtwei = 1+0*p.sqrtwei; % row scaling: make vectors are plain values on bdry
rhs = [p.fillrighthandside; zeros(4*M,1)];   % whose system RHS
kvec = om*exp(1i*p.incang); a = exp(1i*real(conj(kvec) * d));  % Bloch
%s.qpblocha = a;   % for Alpert quad. Not needed w/ Kress reparam

op = [];op.nei = nei;op.buf = op; % dummy qpscatt object for Woods fix one side
pr = qpscatt(du, [], d, op); % omits dd domain
pr.setoverallwavenumber(om); pr.setincidentwave(p.incang); pr.setupbasisdofs;
woods = numel(pr.t.bas)>2;
tu = pr.t; tu.setbloch(a); % use upper qpscatt's qpstrip
for i=1:2, tu.bas{i}.requadrature(M); end; tu.setupbasisdofs;
td = qpstrip(d, omd); %du.isair = 1;
td.addqpftylayerpots(struct('M',M,'nearsing',min(5,omd/2))); td.setbloch(a);
layerextdoms = [du dd]; layerqpstrips = [tu td];
% figure;  layerextdoms.plot; tu.plot; td.plot;

% by hand build system matrix... incl possible Wood's row & col for upper only
tic; Q = [tu.evalbasesdiscrep() zeros(2*M); zeros(2*M,tu.Nf) td.evalbasesdiscrep()];
%nlayers=2; Q = zeros(nlayers * 2*M); % # layers * QP dofs per layer
%for r=1:nlayers, ii = (1:2*M) + (r-1)*2*M;
%  Q(ii,ii) = layerqpstrips(r).evalbasesdiscrep(); end
A = p.fillbcmatrix;    % obstdofs->mismatch block:   self-interaction matrix
for n=[-nei:-1 1:nei]   % direct neighbors
  Ac = p.fillbcmatrix(struct('trans',-d*n));  % don't restrict the doms opt
  A = A + a^n * Ac;
end                    %co = A\rhs(1:p.N); co(1:4), figure; plot(real(co))
% Note: NW and SE blocks of A are I (and -I) if on flat interface (h=0)
% fill B: effect of FTyLPs on each layerextdom-touching segs...
Nq = 2*M; B = [];  % loop over qp-regions, ie layers (UHP, LHP): then segs in p
for r=1:2, Br = zeros(p.N, Nq);  % Nq = total dofs for each regions FTyLP bases
  m=0; Np = numel(vertcat(p.segs.x)); % counter, # pts on all segs in prob
  for i=1:numel(p.segs), si = p.segs(i); ii = m+(1:2*size(si.x,1)); m=ii(end);
    dom = vertcat(si.dom{:}); side = find(vertcat(dom.isair)==r); % side of s
    sa = si.a(side); sb = si.b(side);    % affects layer-touching side only
    t = layerqpstrips(r); nb = numel(t.bas);
    for i=1:nb, b = t.bas{i}; ns = t.basnoff(i)+(1:b.Nf); % air QP bases
      [Bb Bbn] = b.eval(si); %val & nderiv both always need; qpftylp HAZ LR sum!
      Br(ii,ns) = repmat(p.sqrtwei(ii).', [1 b.Nf]) .* [sa * Bb; sb * Bbn];
    end
  end
  B = [B Br];   % stack across for each layer
end
C = [];
for r=1:2, Cr = zeros(2*M, p.N); t = layerqpstrips(r);
  f = utils.copy(t.bas{1}.lp); % dummy qp bas, for transl target L, R evals
  oo = []; oo.dom = layerextdoms(r); 
  for i=1:numel(p.bas)           % all bases in problem or basis obj
    b = p.bas{i}; ns = p.basnoff(i)+(1:b.Nf); % dof indices for this bas
    if utils.isin(layerextdoms(r), b.doms) % restrict to bases affecting airdoms
      %fprintf('i=%d (om=%g), r=%d\n', i, f.om, r)
      for n=nei-2*buf:nei   % general buf case (nei>=buf)
        f.orig = t.Lo-n*d; [Cj Cxj] = b.evalfty(f, oo);
        Cr(:,ns) = Cr(:,ns) + a^n*[Cj;Cxj];
        f.orig = t.Lo-(-1-2*buf-n)*d; [Cj Cxj] = b.evalfty(f, oo);
        Cr(:,ns) = Cr(:,ns) - a^(-1-2*buf-n)*[Cj;Cxj];
      end  % careful, since L's phase = 1 always
    end
  end  % problem basis loop
  C = [C; Cr];    % stack up for each layer
end
%B = 0*B; C = 0*C; % kill qp part for debug
E = [A B; C Q];  % system matrix (no Wood's extra rows)
if woods, [RO RI] = pr.fillbraggamplrow(pr.ikpfix);   % get Bragg ampl matrix
  RI = RI.*repmat(1./sqrt(sum(RI.^2, 2)), [1 size(RI,2)]); % row-normalize
  E = [E(1:p.N+2*M,:); RI zeros(1,td.Nf); E(p.N+2*M+1:end,:)];
  rhs = [rhs; zeros(size(RI,1), 1)];  % insert zero at end since discrep rhs=0
end    % append block row to E and rhs
if v>1,r=sqrt(sum(abs(E).^2, 2));fprintf('row dyn range = %.3g\n',max(r)/min(r))
r=sqrt(sum(abs(E).^2, 1)); fprintf('col dyn range = %.3g\n', max(r)/min(r)),end
if resc, r = sqrt(sum(abs(E).^2, 2)); E = E.*repmat(1./r,[1 size(E,2)]); % rowsc
  rhs = rhs./r; % rhs row scale
  r = sqrt(sum(abs(E).^2, 1)); E = E.*repmat(1./r,[size(E,1) 1]); % col-scale
  %fprintf('cond(E) after rescaling = %.3g\n', cond(E))
end
if 0, % set diag all positive by row rescaling... Fails for multiple segs
  ss=sign(real(diag(E))); E = E.*repmat(ss, [1 size(E,2)]); rhs = rhs.*ss;
end
%figure; plot(eig(E), '+'); axis equal
if v>1, figure; imagesc(real(E)); caxis(0.1*[-1 1]); end % show sys matrix
co = E\rhs;  % solve
fprintf('%dx%d system filled & solved in %.3g s\n',size(E,1),size(E,2),toc);
%co(p.N+1:end) = 0; % kill qp part
if resc, co = co./r.'; end      % undo col scaling

if 1 % check against fresnel flat interface analytic, P(0)
  ko = sqrt(om^2-real(kvec)^2); kod = sqrt(omd^2-real(kvec)^2); sk=ko+kod;
  amplr = (ko-kod)/sk; amplt = 2*ko/sk;
  fprintf('Fresnel flat fracs:  refl, trans = %.15g, %.15g\n\n',amplr.^2,amplt.^2*kod/ko)
end
tic;op = [];op.nei = nei;op.buf = op; % hack the plotting of fields & flux:
pr = qpscatt(du, [], d, op); % dummy qpscatt object for evals: omits dd domain
pr.setoverallwavenumber(om); pr.setincidentwave(p.incang);
pr.t = tu; pr.setupbasisdofs; pr.co = co([1:p.N p.N+[1:tu.Nf]]); % obst + UHP
%pr.pointsolution(pointset(0+0.5i))  % UHP scatt field (zero if omd=om)
du.cloc = nan; dd.cloc = nan;  % excludes strips for plotting
if v, [uu gx gy diu] = pr.showfullfield(struct('ymax', 2.0, 'nx',50));
  title('up field'); close(gcf); end
[bo to bi tiu] = pr.braggampl();
[fuu fdu nu] = pr.braggpowerfracs(struct('table',0)); % heed only fuu
fprintf('Reflected:\tBragg order\tflux fraction\t\ttot: %.14g\n',sum(fuu))
disp([nu fuu])

pr = qpscatt(dd, [], d, op);  % again omit du since confuses the di=d.inside
pr.setoverallwavenumber(om); pr.k = omd; pr.setincidentwave(-acos(real(exp(1i*p.incang))*om/omd)); pr.a = a; % ensure same alpha for LHP
pr.t = td; pr.setupbasisdofs; 
du.k = om; dd.k = om*refr; tu.k = om; td.k = om*refr;
pr.co = co([1:p.N p.N+tu.Nf+[1:td.Nf]]);
%z=-1i; pr.pointsolution(pointset(z))%-p.ui(z)% LHP scatt field: u should be ui
[bo to bid ti] = pr.braggampl();
[fud fdd nd] = pr.braggpowerfracs(struct('table',0,'noinc',1)); % heed fdd
fprintf('Transmitted:\tBragg order\tflux fraction\t\ttot: %.14g\n',sum(fdd))
disp([nd fdd])
fprintf('Err estimate: 2-layer tot flux err = %.3g\n',sum([fuu;fdd])-1)
fprintf('(worst whole-system rad cond mode error = %.3g)\n',max(abs([bid;tiu])))
if v, [ud gx gy did]=pr.showfullfield(struct('ymax', 2.0,'nx',50,'noinc', 1));
  title('down field'); close(gcf); % Now do Fresnel solution:
  if fresnel
    [xx yy] = meshgrid(gx,gy); zz = xx+1i*yy; uex=0*zz; ii = ~isnan(diu);
    uex(ii) = p.ui(zz(ii))+amplr*exp(1i*real(kvec.*zz(ii))); ii = ~isnan(did);
    uex(ii) = amplt*exp(1i*(-kod*imag(zz(ii))+real(kvec)*real(zz(ii))));
  end
  u = zeros(size(uu)); ii = ~isnan(diu); u(ii) = uu(ii); ii = ~isnan(did); u(ii) = ud(ii); % combine the non-nan parts of fields
  up = real(u);
  if fresnel, up = imag(u-uex); end % test against Fresnel soln (vt=0 only!)
  figure; imagesc(gx, gy, up); set(gca,'ydir', 'normal'); colorbar;
  axis equal tight; title(sprintf('u=E_z (TM), \\omega=%.3g n=%.3g d=%.3g \\theta_i=%.3g deg: Reflcoeff=%.3g',om,refr,d,90+p.incang*180/pi,sum(fuu)));
  utils.goodcaxis(up); xlabel('x'); ylabel('y');
  hold on; x = vertcat(s.x); plot([x+d;x;x-d], 'k-');
  fprintf('evaluation and plot time: %.3g s\n',toc)
end
