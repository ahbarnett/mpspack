% Test passing layerpot eval through Jfilter (option of layerpot.eval)
% Seeing if smooth density function can be done with close eval via J filter
% Barnett 1/17/09

clear
verb = 1;            % 0 = text only, 1 = couple figs, 2 = lots of figs
dens = 'd';          % 's' SLP, 'd' DLP
acc = 'm';           % method for acc field eval: c=close-interp, m=many quadr
k = 5;                                          % wavenumber
hlen = 1.0;                                     % half-length of src LP line
y = -1.05; s = segment(70, [-hlen+1i*y hlen+1i*y]);  % line of src pts, fixed y
lp = layerpot(s, dens, k);
R = 1.0;                                         % radius
p = segment(52, [0 R, 0, 2*pi], 'p');           % circle of target pts
if verb>1
  figure; p.plot; hold on; s.plot; axis equal; title('src + targ geom'); end
A = lp.eval(p);                  % default = naive eval
if verb>1, figure; imagesc(real(A)); title('A'); end
%M = ceil(5 + k*R/2 + 40*(k/60)^(1/3));   % estimate max order of J 
o.Jfilter.M = 50;
o.Jfilter.rescale_rad = 1.0;      % radius to scale J-exp coeffs to be O(1) at
o.Jfilter.origin = 0.05+0.07i;    % test displaced origin
B = lp.eval(p, o);               % o opts tells it to use J-filter
fprintf('max abs matrix elementwise err (naive - Jfiltered) = %g\n', ...
        max(abs(A(:)-B(:))))
if acc=='c'
  oc.close = 0.25;                  % use slow interp+gkquad close evals to check
  disp('eval close...'); C = lp.eval(p, oc); disp('done.');
  fprintf('max abs matrix elementwise err (close - Jfiltered) = %g\n', ...
          max(abs(C(:)-B(:))))
end
if verb>1, figure; imagesc(real(B)); title('B'); end
if verb>1, figure; imagesc(log10(abs(B-A))); caxis([-16 0]); colorbar;
  xlabel('src dofs'); ylabel('targ dofs'); title('elementwise LP matrix error')
end

f = @(x) 1+0.5*real(x);          % choose a (smooth) density function
unai = A * f(s.x);               % naive LP eval of field at targs p
ufil = B * f(s.x);               % Jfiltered field eval on targs p

sa = utils.copy(s); sa.requadrature(200); % create >> # src quadr for LP
lpa = layerpot(sa, dens, k);             % make a new LP for later 
if acc=='c'
  uacc = C * f(s.x);
elseif acc=='m'
  Aacc = lpa.eval(p);            % expensive (bigger) matrix fill
  uacc = Aacc * f(sa.x);          % accurate field on targs p
end
fprintf('rel L-infty targ field err (naive - true) = %g\n', ...
        max(abs(unai - uacc))/max(abs(uacc)))
fprintf('rel L-infty targ field err (Jfiltered - true) = %g\n', ...
        max(abs(ufil - uacc))/max(abs(uacc)))

x = -2:0.02:2; [xx yy] = meshgrid(x); g = pointset(xx(:)+1i*yy(:));
po.arrow=0; po.normals=0;
c = (.3+(dens=='d'))/sqrt(k); % plot segment opts, color axis
if verb>0
  u = reshape(lp.eval(g) * f(s.x), size(xx));
  uf = reshape(lp.eval(g, o) * f(s.x), size(xx));
if verb>0, showfield(x, x, real(uf), c, 'Re u, J-filtered');
  hold on; s.plot(1,po); end %p.plot(1,po); end
if verb>0, ua = reshape(lpa.eval(g) * f(sa.x), size(xx)); % slower!
  showfield(x, x, real(ua), c, 'Re u, accurate'); hold on; sa.plot(1,po);
  figure; imagesc(x, x, log10(abs(uf-ua))); caxis([-16 0]); colorbar;
  set(gca,'ydir', 'normal'); title('abs J-filtered log10 err');
  hold on; s.plot(1,po); end
if verb>0, showfield(x, x, real(u), c, 'Re u, naive');
  hold on; s.plot(1,po); p.plot(1,po);
  figure; imagesc(x, x, log10(abs(u-ua))); caxis([-16 0]); colorbar;
  set(gca,'ydir', 'normal'); title('abs naive log10 err');
  hold on; s.plot(1,po); end
end
