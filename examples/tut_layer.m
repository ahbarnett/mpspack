% Example codes from MPSpack tutorial, also generates EPS figures for this doc
% SECTION 6: LAYER POTENTIALS

clear all classes; verb = 0;          % if verb>0, generates EPS figures
tref = segment.radialfunc(100, {@(q) 1 + 0.3*cos(3*q), @(q) -0.9*sin(3*q)});
d = domain([], [], tref, -1); d.k = 10;
d.addlayerpot([], 'd');                    % adds DLP to bdry of d
