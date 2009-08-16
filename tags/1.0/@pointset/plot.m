function h = plot(pts)
% PLOT - plot a pointset
%
%  h = PLOT(pts) plots a pointset with normals, on current figure, returning
%  the graphics handle.
%
% See also: POINTSET

% Copyright (C) 2008, 2009, Alex Barnett, Timo Betcke

g = gcf;
figure(g); hold on;
h = plot(real(pts.x), imag(pts.x), '.');
if ~isempty(pts.nx)                        % show normals...
  l = 0.1;                                 % length of normals
  h = [h; plot([pts.x(:).'; (pts.x(:)+l*pts.nx(:)).'], 'k-')];
end