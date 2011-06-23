function s = smoothnonsym(M, a, b, w)
% SMOOTHNONSYM - non-symmetric oscillatory radial function closed segment
%
%  s = smoothnonsym(M, a, b, w) generates a smooth closed radial function
%   segment with M discretization pts, amplitude a, non-symmetry parameter
%   b, and frequency w (which must be integer).

% Copyright (C) 2011, Alex Barnett, Timo Betcke

s = segment.radialfunc(M, {@(q) 1 + a*cos(w*(q+b*cos(q))), ...
                    @(q) -a*sin(w*(q+b*cos(q))).*w.*(1-b*sin(q)), ...
                    @(q) -a*cos(w*(q+b*cos(q))).*w^2.*(1-b*sin(q)).^2 + ...
                    a*sin(w*(q+b*cos(q))).*w.*b.*cos(q)}); % includes curvature
