% dielectric transmission scattering with Rokhlin hypersingular-cancelling
% barnett 1/12/09

clear classes
verb = 1;                                         % verbosity
k = 8;                                            % overall (ext) wavenumber
n = 1.4;                                            % interior refractive index
M = 200; s = segment.smoothstar(M, 0.2, 3);        % smooth closed segment
di = domain(s, 1); di.setrefractiveindex(n);      % interior
de = domain([], [], s, -1);                       % exterior
o.quad = 'k'; o.ord = 10;             % quadrature for LPs
de.addlayerpotbasis(s, 'd', [], o); de.addlayerpotbasis(s, 's', [], o);
di.addlayerpotbasis(s, 'd', [], o); di.addlayerpotbasis(s, 's', [], o);
setmatch(s, 'diel', 'TM');
pr = scattering(de, di);
if verb, figure; di.plot; hold on; de.plot; axis equal; end

pr.setoverallwavenumber(k);
pr.setincidentwave(pi/2 - pi/20);  % if just angle given, it's a plane wave
pr.fillquadwei; pr.setupbasisdofs; pr.fillbcmatrix; pr.fillrighthandside;
% now hack dofs so sigma,tau apply to both int and ext LPs...
N = pr.N/2; pr.A = pr.A(:,1:N) + pr.A(:,N+1:end);
pr.linsolve; pr.bcresidualnorm, norm(pr.co)
pr.co = [pr.co; pr.co];           % duplicate basis dofs for plotting
pr.pointsolution(pointset(1+1i))        % check u_scatt at one ext pt
if verb,
  opts.dx = 0.05; opts.bb = [-3 3 -3 3]; figure;
  tic; pr.showthreefields(opts); fprintf('\tgrid eval in %.2g sec\n', toc);
end
