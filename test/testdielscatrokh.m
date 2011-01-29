% Dielectric transmission scattering via Rokhlin hypersingular-cancelling scheme
% Comparing Kress hypersingular spectral vs Kapur-Rokhlin & Alpert quadratures.
% 1/12/09, updated for many doms per bas object, 6/2/10, Alpert 1/29/11
% Alex Barnett
clear all classes; verb = 1;                      % verbosity (0=text, 1=plot)
k = 8;                                            % overall (ext) wavenumber
n = 1.4;                                          % interior refractive index
M = 110; s = segment.smoothstar(M, 0.2, 3);       % smooth closed segment
di = domain(s, 1); di.setrefractiveindex(n);      % interior
de = domain([], [], s, -1);                       % exterior
o.quad = 'a'; o.ord=16;                           % quadr scheme and order
s.addinoutlayerpots('d', o);                      % new double-sided layerpot
s.addinoutlayerpots('s', o);                      % "
setmatch(s, 'diel', 'TM');
pr = scattering(de, di);
if verb, figure; di.plot; hold on; de.plot; axis equal; end

pr.setoverallwavenumber(k);
pr.setincidentwave(pi/2 - pi/20);  % if just angle given, it's a plane wave
pr.fillquadwei; pr.setupbasisdofs;
pr.fillrighthandside;
pr.fillbcmatrix;
pr.linsolve;
fprintf('resid = %.3g, coeff nrm = %.3g\n', pr.bcresidualnorm, norm(pr.co))
u = pr.pointsolution(pointset(1+1i))        % check u_scatt at one ext pt
ugood =  1.176452635715030 - 0.798366817843056i % (Kress M=110 to 1e-15)
fprintf('error at one point = %.3g\n', abs(u-ugood))
% Note Kapur-Rokh needs N=600 to get ans to 1e-8, N=400 to 1e-6. Terrible.
% Alpert 16th gets 1e-11 at N=110 (where Kress is 1e-15). Not bad.
if verb, opts.dx = 0.05; opts.bb = [-3 3 -3 3]; figure;
  tic; pr.showthreefields(opts); fprintf('\tgrid eval in %.2g sec\n', toc);
end
