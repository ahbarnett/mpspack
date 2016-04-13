% Inclusion demo for ninc paper, cleaned up from fig_tsweeps. Barnett 12/6/15
% edited for mpspack, 4/13/16.
%
% If kex ("exact" wavenumber) set, skips the minimization search of tension.
% If swp=1, just does sweep and saves data used by fig_tsweeps.m
clear; verb = 1;  % verbosity; 0 = no mode, 1 = plot the mode

krng = 1;  % choose k range to explore (1=medium, 1st row of Table; 2=high)
swp = 0;  % 0: no sweep (do inclusion bound), 1: write sweeps for fig_tsweeps.m

if krng==1
  kwin = [40.50 40.55];  % since another at 40.5679
  %kex = 40.51282199500848; % use to override search; exact ev, eg from fred-det
  kswp = 40:0.001:41;      % 1e3 steps; Erng = max(kswp)^2-min(kswp)^2 = 81
  %kswp = 40.50:0.005:40.55;  % to look at
  M = 700;      % bdry
  N = 350;      % basis
  tau = 0.025; %0.016;  % MFS dist
  dx = 0.01; cax = 2;   % plot grid, color scale
else
  kwin = [405.000 405.005];
  kex = 405.0032695182281;        % 3e-16 error if change last digit by 1
  kswp = 405:0.0001:405.1;    % 1e3 steps, same E range as above.
  M = 5000;      % bdry
  N = 2500;      % basis
  tau = 0.004;
  if swp, M = 2500; N = 1500; tau=0.005; end %  takes 5 hrs for 1e3 steps.
  dx = 0.002; cax = 2;   % plot grid, color scale
end

% geom setup
s = segment.smoothnonsym(M, 0.3, 0.2, 3);  % shape
d = domain(s,1);
L = d.perim; ktyp = max(kwin);
fprintf('mean M ppw = %.3g; N ppw = %.3g\n',M/(ktyp*L/2/pi), N/(ktyp*L/2/pi))
dE = 4*pi/d.area;                   % Weyl mean level spacing in E=k^2
fprintf('mean dE spacing %.3g, dk spacing %.3g\n',dE, dE/2/ktyp)

% basis setup
d.addmfsbasis(s, [], struct('tau',-tau)); d.bas{1}.realflag=1; % real for speed
p = evp(d);       % create eigenvalue problem
p.updateN(N); if verb>2, figure; d.plot; d.showbasesgeom; end
% method used for computing t for search, and getting best coeff vec y:
%to = struct('ten','h', 'reg','t', 'eps',1e-14, 'Fh',0);  % classical ttilde
%to = struct('ten','h', 'reg','t', 'eps',1e-14, 'Fh',1);  % our new ttilde
%to = struct('ten','b', 'reg','t', 'eps',1e-14, 'Geps',1e-12, 'Fh',1);  % our new ttilde; note Geps can't be too small
to = struct('ten','b', 'reg','s', 'eps',1e-14, 'Geps',1e-12, 'Fh',1);  % our new ttilde; note Geps can't be too small
%to = struct('ten','v', 'reg','t', 'eps',1e-14, 'Fh',0);
% Observe that 's' is no slower than 't', so stick to the plain SVD regularize.


if swp  % .................. sweep only .............................
  ns = 100; ts = nan(ns,numel(kswp));  % tension data
  for i=1:numel(kswp), k=kswp(i);
    tic; t = evp.gsvdtension(d, k, to);
    nkeep = min(ns,numel(t)); ts(1:nkeep,i) = t(1:nkeep);
    fprintf('k=%.15g: \tt_min = %.3g   (%.3g sec)\n',k,min(t),toc)
  end
  figure; plot(kswp.^2,min(ts,[],1),'-');
    hold on; plot(kswp.^2,ts,'b.'); % all curves
  xlabel('E'); ylabel('$\tilde{t}_h$','interpreter','latex');
  axis([min(kswp)^2 max(kswp)^2 0 4])
  %save rfn_newtsweep_40k41       % read by fig_tsweeps.m
  %save rfn_newtsweep_405k405.1
  stop
end

% if kex not set, use recursive search on min gen sing val to find best k...
%profile clear; profile on
if ~exist('kex','var')   % sing val minimization search, equal grid inside kwin
  fprintf('finding plain val-part Neu norm tension minima in [%.16g,%.16g]...\n',kwin(1),kwin(2))
  kgrid = linspace(kwin(1),kwin(2),5);     % 3 is too small, can jump away
  io = []; io.xtol = 1e-12; io.maxslope = ktyp/2; io.verb = 1; % gridminfit opts
  nsingvalskeep = 10;
  tfunc = @(k) utils.lowestn(evp.gsvdtension(d, k, to),nsingvalskeep);
  tic, [kj tj mininfo] = evp.gridminfit(tfunc, kgrid, io);   % the meat
  fprintf('min search (%d evals) in %g s\n',numel(mininfo.xs),toc)
  if verb, figure; plot(mininfo.xs,mininfo.ys,'+'); hold on;
    plot(kj,min(tj,[],1),'r*'); title('gridminfit: tension evals vs k'); end
  for j=1:numel(kj), fprintf('k_%d = %.16g, \tt = %.3g\n',j,kj(j),min(tj(:,j))), end
  j = 1; kj = kj(j); tj = min(tj(:,j)); fprintf('keeping j=%d\n',j)
else
  kj = kex; fprintf('using k_j = k_ex = %.15g\n',kj);    % override w/ known k
end
%profile viewer

d.k = kj;  % from this point, kj is the only info used from the minimization
E = kj^2; fprintf('Weyl j est = %d\n',round(E*d.area/(4*pi) - L/(4*pi)*sqrt(E)))

% Now the only piece of info going forward is kj, the estimated eigenwavenumber
% Compute best coeff vec y at k_j:
tic, [ts,V] = evp.gsvdtension(d, kj, to); fprintf('tens+coeffs in %.3g s\n',toc)
l = find(abs(ts)==min(abs(ts))); t = ts(l); y = V(:,l);
fprintf('coeff vec ||y|| = %.3g\n',norm(y))
fprintf('t from tension used for search = %.3g\n',t)

[Ab Anb] = d.bas{1}.eval(s); % get self-eval matrices at k_j, O(N^2) but acc
unb = Anb*y; ub = Ab*y;      % eval u,u_n on bdry (cf FMM loses digits, 1e-12)
D = 2*pi*repmat(1./d.speed,[1 M]).*circulant(quadr.perispecdiffrow(M)); % d/dt
utb = D*ub;                  % u_t on bdry

xn = real(conj(s.x).*s.nx);  % check interior norm (assuming Neu BC)
integrand = (1/2)*xn.*(abs(ub).^2 - (abs(utb)).^2/kj^2); %int nrm form, Neu case
intnrm = sqrt(sum(s.w'.*integrand));
fprintf('\tinterior norm of u from coeff vec y = %.16g\n',intnrm)
y = y/intnrm; ub = ub/intnrm; unb = unb/intnrm; utb = utb/intnrm;  % normalize

ten = sqrt(sum(s.w'.*abs(unb).^2)); fprintf('classical ttilde = %.3g\n',ten)

if verb>1, % make mode figure
  fprintf('grid eval at dx=%g (ppw=%.3g)...\n',dx,2*pi/(ktyp*dx))
  p.co = y; figure; tic
  [u gx gy di] = p.showsolution(struct('dx', dx, 'fmm', 1)); % fast
  fprintf('grid eval in %.3g s\n',toc), caxis(cax*[-1 1]);
  ii = find(di==1);    % interior indices in u array
  unrmdx = dx*norm(u(ii)); fprintf('check crude int nrm = %g\n',unrmdx)
  box off, axis off, title '', colorbar off
  %print -depsc2 neu_mode_high.eps   % resolution 0.002 could be improved
end

Ej = kj^2; Cen = Cennenbach(p); Eerrclas = Cen*Ej*ten % classical incl bnd on E
fprintf('Eerrclas / E = %.3g\n', Eerrclas/Ej)

Fh = evp.spectralfiltermatrix(s,kj);   % Our bounds w/ F_h(..) filter applied:
Funb = Fh*unb;                       % filtered bdry u_n
ttilde = sqrt(sum(s.w'.*abs(Funb).^2))
fprintf('rel err ttilde/E = %.3g\n',ttilde/E)
Cest = 1.6;
fprintf('rel bnd Cest.ttilde/E = %.3g\n',Cest*ttilde/E)
if verb>1, figure; plot(cumsum(s.w),[unb Funb],'.-'); title('u_n and F_h u_n err funcs'); end

% predict nearby E_j:
%typslope = 0.65; kslope = typslope*2*kj
%predkdists = ts(1:10)/kslope
