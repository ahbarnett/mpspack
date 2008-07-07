% PLOT - plot a pointset
%
%  h = PLOT(pts) plots a pointset with normals, on current figure, returning
%   the graphics handle.
%
% See also: POINTSET

function h = plot(pts)

g = gcf;
figure(g); hold on;
h = plot(real(pts.x), imag(pts.x), '.');
% show normals...
l = 0.1;                                 % length of normals
h = [h; plot([pts.x(:).'; (pts.x(:)+l*pts.nx(:)).'], 'k-')];
