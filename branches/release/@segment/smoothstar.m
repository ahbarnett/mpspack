function s = smoothstar(M, a, w)
% SMOOTHSTAR - single-freq oscillatory radial function closed segment
%
%  s = smoothstar(M, a, w) generates a smooth closed radial function segment
%   with M discretization pts, amplitude a, and frequency w (integer)
%
% Copyright (C) 2008, 2009, Timo Betcke, Alex Barnett

s = segment.radialfunc(M, {@(t) 1 + a*cos(w*t), @(t) -w*a*sin(w*t), ...
                    @(t) -w^2*a*cos(w*t)});
