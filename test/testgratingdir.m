% basic Dirichlet diffraction grating, connected surface, unwrapped Alpert quad
% barnett 1/29/11

clear; v = 1;        % verbosity: 0 no pics, 1 final pic, 2 diagnostics
d = 1.0;   % problem x-periodicity 
vt = -0.2;  % vertical offset
h=.4; s = segment(50, {@(t) 1i*vt + .5-t+1i*h*sin(2*pi*t), @(t) -1+2i*pi*h*cos(2*pi*t), @(t) -4i*pi^2*h*sin(2*pi*t)},'p'); % wobbly wall (note R-to-L to get normal up)
om = 10;                                        % incident wavenumber
o.nei = 1; o.buf = 0; o.M = 90;                 % M = # nodes per FTy LP

de = domain([], [], s, -1); de.cloc=0; % hack so (most of) R2 in this domain,
% in particular the B segment is inside, so fillbraggamplrow can work.
%de = qpstrip(d, om, struct('seg',s, 'pm', -1)); de.exterior = 1; s.dom{1} = de;
o.quad = 'a'; o.ord = 16;       % quadrature
de.addlayerpot(s, [-1i*om 1], o);              % adds CFIE to segment
s.setbc(1, 'D', []);                           % homog Dirichlet BCs
p = qpscatt(de, [], d, o);                     % create problem instance
p.setoverallwavenumber(om);
%p.setincidentwave(-pi/2);                      % normal inc, alpha=1 for debug
p.setincidentwave(-pi/5);                      % generic non-Wood's angle
%p.setincidentwave(-acos(1-2*pi/om)-0*1e-14);   % single Wood's anomaly
if v, figure; p.plot; hold on; p.showbragg; end
s.qpblocha = p.a;   % tell segment its Bloch alpha (put this into setinc...?)
p.fillquadwei;             % this is only for obstacle mismatch (blocks A,B)
p.sqrtwei = 1+0*p.sqrtwei; % row scaling: make vectors are plain values on bdry
p.fillrighthandside;
tic; p.fillbcmatrix; toc
tic; p.co = p.linsolve; toc;    % least-squares soln
z = pointset(0.2+0.5i); tic; p.pointsolution(z), toc
de.cloc=0; [u d n] = p.braggpowerfracs(struct('table',1)); de.cloc = nan; % hack which returns domain.inside to its original behavior for fillbraggamplrow
fprintf('up flux err=%.3g\n',1-sum(u))

c = p.co(1:p.N);                     % eta
if v>1, figure; plot([s.t; 1+s.t], abs([c; c]),  '+-'); % abs(dens) is periodic
title(sprintf('P=%d Z=%d',o.nei,o.buf)); end

de.cloc = nan; % needed to tell de.inside its grating, exclude lower halfplane:
if v, tic; [u gx gy di] = p.showfullfield(struct('ymax', 1.3)); toc, end
