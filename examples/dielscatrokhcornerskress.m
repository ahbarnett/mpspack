% dielectric transmission scattering with Rokhlin hypersingular-cancel scheme
% and trying Kress-reparametrized corners. Various shapes: tri, cuspy astroid...
% Barnett 3/16/13

clear; verb = 1;                      % verbosity
k = 10;                                           % overall (ext) wavenumber
n = 1.5;                                  % interior refractive index (real)
Ns = 40:20:180;   % shows 120 enough for kressq=6 @ 1e-11 error
for i=1:numel(Ns), N = Ns(i); % convergence .........
  % pick shape by uncommenting line(s):
  %s = segment.smoothstar(N, 0.3, 3);     % smooth closed segment
  %o.kressq=6; s = segment.polyseglist(N, [1,-1+1i,-.5-1i], 'pc', o); %tri
  %o.kressq=6; s = segment.polyseglist(N, [1,-1+1i,1-0.5i], 'pc', o); %thin tri
  %o.kressq=8; s = segment(N, [1+1i 1 3*pi/2 pi], 'pc', o); s = [s s.rotate(pi/2) s.rotate(pi) s.rotate(3*pi/2)]; % 4-cusp astroid shape
  r=0.001; a = sqrt(r*(2+r)); t=atan(a); % rounding radius of rounded astroid
  o.kressq=5; ss(1) = segment(N/5, [1-a,r,-pi/2+t,pi/2-t], 'pc', o); % rounded
  o.kressq=5; ss(2) = segment(N, [1+1i 1 3*pi/2-t pi+t], 'pc', o); % main curve
  s = [ss ss.rotate(pi/2) ss.rotate(pi) ss.rotate(3*pi/2)]; clear ss
di = domain(s, 1); di.setrefractiveindex(n);      % interior
de = domain([], [], s(end:-1:1), -1);             % exterior (reversed order)
o.quad = 'm';                                     % Kress spectral quadr corr
s.addinoutlayerpots('d', o);                      % double-sided layerpot
s.addinoutlayerpots('s', o);                      % "
setmatch(s, 'diel', 'TM');
pr = scattering(de, di);
if verb & i==1, figure; di.plot; hold on; de.plot; axis equal; drawnow, end

pr.setoverallwavenumber(k);
pr.setincidentwave(pi/6);  % if just angle given, it's a plane wave
pr.fillquadwei; pr.setupbasisdofs;
pr.sqrtwei = 1+0*pr.sqrtwei; % unweight (kill Bremer), is better for GMRES
pr.fillrighthandside;
pr.fillbcmatrix;
%si=sign(real(diag(pr.A))); pr.A = pr.A.*repmat(si, [1 size(pr.A,2)]); pr.rhs = pr.rhs.*si; % make all diag signs the same, makes Id+cpt not +-Id+cpt.
pr.linsolve;
u=pr.pointsolution(pointset(1+1i));        % check u_scatt at one exterior pt
fprintf('N=%d:\tu(pt) = %.16g + %.16gi\n',N,real(u),imag(u))
end   % .............

if verb, opts = []; opts.dx = 0.03; opts.bb = [-2 2 -2 2]; figure;
  tic; pr.showthreefields(opts); fprintf('\tgrid eval in %.2g sec\n', toc);
  hold on; plot([vertcat(s.x);s(1).x(1)],'k-');
end
