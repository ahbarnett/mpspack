% grating of dielectric obstacles via MPSpack to 14 digits. Barnett 6/10/10
d = 1.0;                                            % problem period
N = 80; s = scale(segment.smoothstar(N, 0.3, 3), 0.35); % smooth closed curve
o.nei = 2; o.buf = 1; o.M = 150; % method params, M = # nodes per FTy LPs
di = domain(s, 1); di.setrefractiveindex(1.5);      % obstacle, refractive index
de = domain([], [], s, -1);                         % obstacle's exterior
s.addinoutlayerpots('d'); s.addinoutlayerpots('s'); % add Rokhlin LP scheme 
s.setmatch('diel', 'TM');                           % TM dielectric continuity
p = qpscatt(de, di, d, o);                          % set up the problem
p.setoverallwavenumber(10);                         % incident wavenumber
p.setincidentwave(-pi/5);                           % incident plane wave angle
p.solvecoeffs;                                      % fill matrices and solve
[u d n] = p.braggpowerfracs(struct('table',1));     % intensity of Bragg orders
p.showfullfield(struct('ymax', 1));                 % plot Re part of full field
