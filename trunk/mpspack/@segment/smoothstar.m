% SMOOTHSTAR - generate oscillatory radial function segment
%
%  s = smoothstar(M, a, w) generates a smooth closed radial function segment
%   with M discretization pts, amplitude a, and frequency w (integer)

function s = smoothstar(M, a, w)

R = @(t) 1 + a*cos(w*t); Rt = @(t) -w*a*sin(w*t);
Z = @(s) exp(2i*pi*s).*R(2*pi*s);
s = segment(M, {Z, @(s) 2*pi*(1i*Z(s) + exp(2i*pi*s).*Rt(2*pi*s))}, 'p');
