% Compute Dirichlet eigenvalues of a 40-sided polygonal approximation to the UK
% Dave Hewett 3/11/16, tweaked by Alex Barnett.
% based on tut_evp.m from mpspack examples. Needs hewett_ukmap_points.m

% Note by Barnett: due to nasty corners, this is not a terribly impressive
% demo of MPSpack, but we include it for educational purposes.
% It takes around 10 minutes to compute 7 eigenvalues to around 4-5 digit
% accuracy. Partly this is due to MATLAB overhead in the 40^2 segment
% interactions; partly due to needing a large number of points per segment
% due to non-geometric corner quadratures.
clear all classes; verb = 0;
o.kressq = 6;    % corner-packing parameter for Kress reparametrization
uk_pts = flipud(hewett_ukmap_points(40)).';
n = 60;
s = segment.polyseglist(n, uk_pts, 'pc', o);  % UK
d = domain(s, 1);        % create an interior domain
s.setbc(-1, 'D');        % Dirichlet BC's applied on inside of segment
p = evp(d);              % sets up eigenvalue problem object
d.addlayerpot(s, 'd');          % DLP basis set appropriate for Dir BC

%profile clear; profile on
o.tol = 1e-2; o.modes = 1; tic; p.solvespectrum([5 15], 'fd', o); toc
%profile off; profile viewer

o = []; o.dx = 0.01; o.FMM=1; tic; p.showmodes(o); toc
%print -dpng ../gallery/hewett_uk_modes.png
