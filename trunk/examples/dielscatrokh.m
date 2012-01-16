% dielectric transmission scattering with Rokhlin hypersingular-cancel scheme
% Barnett 6/18/10

clear all classes; verb = 1;                      % verbosity
k = 30;                                           % overall (ext) wavenumber
n = 1.5;                                          % interior refractive index
M = 330;                                          % gives 13 digit acc at k=30
s = segment.smoothstar(M, 0.3, 3);                % smooth closed segment
di = domain(s, 1); di.setrefractiveindex(n);      % interior
de = domain([], [], s, -1);                       % exterior
o.quad = 'm';                                     % Kress spectral quadr
s.addinoutlayerpots('d', o);                      % new double-sided layerpot
s.addinoutlayerpots('s', o);                      % "
setmatch(s, 'diel', 'TM');
pr = scattering(de, di);
if verb, figure; di.plot; hold on; de.plot; axis equal; end

pr.setoverallwavenumber(k);
pr.setincidentwave(pi/6);  % if just angle given, it's a plane wave
pr.fillquadwei; pr.setupbasisdofs;
pr.fillrighthandside;
pr.fillbcmatrix;
pr.linsolve;
pr.pointsolution(pointset(1+1i))        % check u_scatt at one exterior pt
if verb, opts.dx = 0.03; opts.bb = [-2 2 -2 2]; figure;
  tic; pr.showthreefields(opts); fprintf('\tgrid eval in %.2g sec\n', toc);
end
