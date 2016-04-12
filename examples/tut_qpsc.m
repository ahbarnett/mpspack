% Example codes from MPSpack tutorial, also generates EPS figures for this doc
% SECTION 8: PERIODIC SCATTERING.    Barnett 6/20/10
verb = 0;                      % verb>0 generates EPS figures

% Grating of dielectric obstacles, to 14 digits accuracy, low wavenumber
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
p.pointsolution(pointset(0.3+0.7i))                 % eval full field at a point

if verb, set(gcf,'paperposition', [0 0 4 3]); axis off; colorbar off; title('');
  p.showbdry(struct('arrow',0,'normals',0,'label',0,'blobs',0));
  print -depsc2 ../doc/figs/qpsc.eps; end

% More complicated case of dielectric-coated Dirichlet obstacle with
% nearby Dirichlet disc obstacle which wraps around unit cell bdry (13 digits)
clear; d=1; s = scale(segment.smoothstar(130, 0.3, 3), 0.35);   % as above
sm = s.scale(0.6); di = domain(s, 1, sm, -1);       % diel domain
di.setrefractiveindex(1.5);
c = segment(60, [.5-.6i, .2 0 2*pi]);          % small circle segment
de = domain([], [], {s c}, {-1 -1});          % twice-punctured exterior domain
s.addinoutlayerpots('d'); s.addinoutlayerpots('s'); % as above
om = 20;                                            % incident wavenumber
di.addlayerpot(sm, [-1i*om 1]); de.addlayerpot(c, [-1i*om 1]); % CFIEs for Dir
s.setmatch('diel', 'TM'); sm.setbc(1, 'D', []); c.setbc(1, 'D', []); % BCs
o.nei = 2; o.buf = 1; o.M = 150; p = qpscatt(de, di, d, o);
p.setoverallwavenumber(om); p.setincidentwave(-acos(1-2*pi/om));
p.solvecoeffs;                                      % fill matrices and solve
[u d n] = p.braggpowerfracs(struct('table',1));     % intensity of Bragg orders
p.showfullfield(struct('ymax', 1.2));

if verb, axis off; colorbar off; title('');
  p.showbdry(struct('arrow',0,'normals',0,'label',0,'blobs',0));
  set(gcf,'paperposition', [0 0 4 3]);
  print -depsc2 ../doc/figs/qpsc_coated.eps; end
