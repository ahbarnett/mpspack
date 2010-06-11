% use Dirichlet and diel. Rokh scatt as examples. Also tests domainindices
% barnett 6/10/10

clear all classes; sys = 't'; % system 'd' or 't' for scattering
k = 8;                                            % overall (ext) wavenumber
n = 1.4;                                          % interior refractive index
M = 100; s = segment.smoothstar(M, 0.2, 3);       % smooth closed segment
de = domain([], [], s, -1);                       % exterior
if sys=='d'
  de.addlayerpot(s, [-1i*k 1]); s.setbc(1, 'D', []); % Dir BCs and CFIE form
  pr = scattering(de, []);
elseif sys=='t'
  di = domain(s, 1); di.setrefractiveindex(n);    % interior
  s.addinoutlayerpots('d');                       % new double-sided layerpot
  s.addinoutlayerpots('s');                       % "
  setmatch(s, 'diel', 'TM'); pr = scattering(de, di);
end
pr.setoverallwavenumber(k);
pr.setincidentwave(pi/2 - pi/20);  % if just angle given, it's a plane wave
pr.solvecoeffs;
p = pointset([1+1i;0+.3i]); pr.pointsolution(p)   % u_scatt at ext & int pts
A = pr.evalbases(p); u = A*pr.co         % note get zero rather than NaN
disp('u should agree with the direct pointsolution case');
p.nx = [1;1]; [A Ax Ay] = pr.evalbases(p);     % check handles normals too
pr.domainindices(p) % domain indices should be 1,2 (ext, int)

% test the opts.dom override (for sys='t')...
[A Ax Ay] =  pr.evalbases(p, struct('dom', de)); u = A*pr.co

