%  PLOT - shows a segment
%
% h = PLOT(seg) plots a segment, quadrature points, etc.
%
% h = PLOT(seg, pm) plots segment if pm=1, or its reversal is pm=-1

function h = plot(s, pm)

if nargin<2, pm = 1; end

g = gcf;
figure(g); hold on;
h = plot(real(s.x), imag(s.x), '.-');
% show normals
l = 0.1;                                 % length of normals
plot([s.x(:).'; (s.x(:)+l*pm*s.nx(:)).'], 'k-');   % uses sign from pm
