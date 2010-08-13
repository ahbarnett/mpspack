% basic test routine for EVP class in MPSpack. Barnett 8/13/10

clear all classes; a = 0.3; b = 0.2; w = 3;  % shape params, smooth closed curve
s = segment.radialfunc(160, {@(q) 1 + a*cos(w*(q+b*cos(q))), ...
                    @(q) -a*sin(w*(q+b*cos(q))).*w.*(1-b*sin(q)), ...
                    @(q) -a*cos(w*(q+b*cos(q))).*w^2.*(1-b*sin(q)).^2 + ...
                    a*sin(w*(q+b*cos(q))).*w.*b.*cos(q)});  % includes curvature
d = domain(s,1);
s.setbc(-1, 'D');             % homog dirichlet BCs (on inside: note -1)
d.addlayerpot(s, 'd');
p = evp(d);
p.solvespectrum([2 10]);
disp('eigenwavenumbers (sqrt of eigenvalues):'), p.kj
