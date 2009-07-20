% Example codes from MPSpack tutorial, also generates EPS figures for this doc
% SECTION 4: EXTERIOR MFS

cd ../doc                             % so figures write out to doc/ directory
clear all classes; verb = 0;          % if verb>0, generates EPS figures

% exterior domain
tref = segment.radialfunc(50, {@(q) 1 + 0.3*cos(3*q), @(q) -0.9*sin(3*q)});
d = domain([], [], tref, -1);
   
% TIMO add code for mfs

%...

% multiply-connected domains. 1 hole...
tref.disconnect;                         % clears any domains from segment
c = segment([], [0.5 0.4 0 2*pi]);       % new circular segment
d = domain(tref, 1, c, -1);
% 2 holes...
tref.disconnect; c.disconnect;
smtref = tref.scale(0.3);                % create new rescaled copy of tref
smtref.translate(-0.3+0.5i);             % move the segment smtref
d = domain(tref, 1, {c smtref}, {-1 -1});
if verb  % generate f:doms a
  figure; set(gca, 'fontsize', 14); d.plot; axis off;
  print -depsc2 twoholes.eps
end
  
