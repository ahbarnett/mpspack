function s = smoothstar(M, a, w, p)
% SMOOTHSTAR - single-freq oscillatory radial function closed segment
%
%  s = smoothstar(M, a, w) generates a smooth closed radial function segment
%   with M discretization pts, amplitude a, and frequency w (integer)
%
%  s = smoothstar(M, a, w, p) rotates it by an angle offset of p radians.

% Copyright (C) 2008 - 2012, Alex Barnett, Timo Betcke

if nargin<4, p = 0; end     % default angle offset
s = segment.radialfunc(M, {@(t) 1 + a*cos(w*(t-p)), @(t) -w*a*sin(w*(t-p)), ...
                    @(t) -w^2*a*cos(w*(t-p))});
