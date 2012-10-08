% Example codes from MPSpack tutorial, also generates EPS figures for this doc
% FMM ACCELERATED LAYER POTENTIALS, for exterior Dirichlet demo. Barnett 3/12/11

clear; s = segment.smoothstar(400, 0.3, 9);
d = domain([], [], s, -1); d.k = 10;
d.addlayerpot(s, 'D');                     % adds DLP to segment s
f = @(z) besselh(0,d.k * abs(z-0.3-0.2i)); % known exterior field
s.setbc(1, 'D', [], @(t) f(s.Z(t)));       % its Dirichlet data
p = bvp(d);


p.solvecoeffs;

z = 1+1.2i; err = abs(f(z) - p.pointsolution(pointset(z)))
