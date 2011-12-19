% demo point source incident wave scattering from two dielectric discs.
% For Hakan Tureci. Barnett 12/17/11
clear all classes;

M=70; R=0.8; n=2.0; k=9.8; % R=radius. 4M total unknowns
s = [segment(M, [-1 R 0 2*pi]), segment(M, [1 R 0 2*pi])]; % discs at +-1
d(1) = domain(s(1), 1); d(1).setrefractiveindex(n);
d(2) = domain(s(2), 1); d(2).setrefractiveindex(n);
e = domain([], [], {s(1) s(2)}, {-1 -1});
s.setmatch('diel','tm');
o.quad = 'm';   % choose Kress spectral quadrature
s.addinoutlayerpots('d', o); s.addinoutlayerpots('s', o);
p = scattering(e,d);
p.setoverallwavenumber(k);
x0 = 0; p.setincidentwave(x0, 'pt'); % point src @ x0
%p.setincidentwave(pi/3);            % plane wave
tic; p.solvecoeffs; toc
p.pointsolution(pointset(0.2+1i))  % eval scatt field at a pt (M=70: 1e-7 err)
tic; o.bb = [-1 5 -1 3]; o.dx=0.03; o.sepfigs=1; p.showthreefields(o); toc
