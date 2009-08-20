% dielectric transmission scattering with Rokhlin hypersingular-cancelling
% Comparing Kress hypersingular spectral (great) vs Kapur-Rokhlin quadrature
% Barnett 1/12/09

clear all classes
verb = 1;                                         % verbosity
k = 8;                                            % overall (ext) wavenumber
n = 1.4;                                          % interior refractive index
M = 110; s = segment.smoothstar(M, 0.2, 3);        % smooth closed segment
di = domain(s, 1); di.setrefractiveindex(n);      % interior
de = domain([], [], s, -1);                       % exterior
o.quad = 'm';                          % Kress spectral quadr
%o.quad = 'k'; o.ord = 10;             % Kapur-Rokh quadrature for LPs, M=200 ok
de.addlayerpot(s, 'd', o); de.addlayerpot(s, 's', o);
di.addlayerpot(s, 'd', o); di.addlayerpot(s, 's', o);
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
% compare: 1.176452635715030 - 0.798366817843056i   (Kress M=110 to 1e-15)
% Note Kapur-Rokh needs N=600 to get ans to 1e-8, N=400 to 1e-6. Terrible.
if verb, opts.dx = 0.05; opts.bb = [-3 3 -3 3]; figure;
  tic; pr.showthreefields(opts); fprintf('\tgrid eval in %.2g sec\n', toc);
end
