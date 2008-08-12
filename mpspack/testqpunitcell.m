% test photonic bands QP unit cell object, basis discrep evaluation, QP solve
% barnett 8/6/08

clear all classes
qpb = 'j';                         % basis for doing QP: j, pw, l, lc, ...
verb = 0;                          % 0=no figs, 1=final fig; 2=also geom figs
uc = qpunitcell(1, .5+1i, [], 20);    % M controls # LP quad pts too     .5+1i
if verb>1, uc.setbloch(1, .5); end                % should give warning
uc.setbloch(-1, -1);               % furthest point from origin in Brillouin
if verb>1, uc.plot; figure; uc.showbrillouin; end

k = 10;
%N = 10 + ceil(k*uc.diam/2 + 40*(k/60)^(1/3)); % # dofs/2, from latticesum code
N = 30;
if qpb=='j'
  %opts.rescale_rad = .7; uc.addregfbbasis(0, N, k, opts);
  uc.addregfbbasis(0, N, k);
elseif qpb=='pw'
  uc.addrpwbasis(N, k);
elseif qpb=='l'      % ........ L+B+R+T SLP alone
  uc.addlayerpotbasis(uc.seg, 's', k);   % note SLP alone beats SLP + DLP!
%  uc.addlayerpotbasis(uc.seg, 'd', k);  % bad since hypersingular T self-int
elseif qpb=='lc'     % ........ L+B + sticking-out density copies.
  uc.addqpuclayerpots(k);
end
uc.setbloch(1+2i); %(1+1i); % or use uc.setbloch(0) for test Im total = 0
tic; Q = uc.evalbasesdiscrep; toc
%profile clear; profile on; Q = uc.evalbasesdiscrep; profile off; profile viewer
fprintf('cond(Q) = %.2g\n', cond(Q))
b = mfsbasis(1+0.1i, [], [], k);       % make ext source as a single MFS charge
%b = mfsbasis(0.2+0.1i, [], [], k);     % make int source as a single MFS charge
di = b.evalunitcelldiscrep(uc);        % discrep col vec
tic; co = -Q\di; toc                   % uc basis coeffs to cancel discrep
%tic; co = gmres(Q,di,[],[],30); toc    % GMRES w/ maxit
pr = bvp(uc); pr.co = co;              % set up dummy bvp, pass in co, to plot
fprintf('coeff nrm = %.2g, residual l2 error = %.2g\n',norm(co),norm(Q*co + di))
o.dx = .02; [u gx gy dii] = pr.gridsolution(o);
[xx yy] = meshgrid(gx, gy); zz = xx+1i*yy;
us = NaN*zeros(size(u)); ii=find(~isnan(dii));
us(ii) = b.eval(pointset(zz(ii)));     % compute src field inside on same grid
if verb, showfield(gx, gy, u + us, [], 'src + QP field'); end
iii = find(uc.inside(2*zz));           % points inside central 1/4 of domain
% This should be small if a reg Helm soln was used as us...
fprintf('L2 interior domain norm u + us = %.2g\n', norm(u(iii)+us(iii))*o.dx)
if verb
  figure; imagesc(log10(abs(Q))); caxis([-16 1]); colormap(jet(256)); colorbar;
end
if verb % compare against 'j' default case...
  uc.clearbases; uc.addregfbbasis(0, 30, k);
  Qr = uc.evalbasesdiscrep; 
  co = -Qr\di; pr.co = co; [uref] = pr.gridsolution(o);
  fprintf('L2 domain error from J calc = %.2g\n', norm(u(ii) - uref(ii))*o.dx)
  if verb, showfield(gx, gy, u - uref, [], 'QP field err from J calc'); end
end
