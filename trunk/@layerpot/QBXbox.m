function [box i] = QBXbox(b, t, h, p, o)
%
% box = QBXbox(b, t, h) returns in box.v vertices of approximating polygon for
%  a QBX close-eval box centered a param t, on the left side of a segment.
%  h>0 is the box scale, in [0,1]-scaled segment parameter.
%  If h has 2 elements, then first is real-half-width, second is imag-height,
%  allowing aspect ratio to change from 2x1.
%
% box = QBXbox(b, t, h, [], o) controls opts ....
%
% [box i] = QBXbox(b, t, h, p, o) also returns logical array i corresponding
%  to whether each element of list of complex numbers p is inside the box.
%
% TO DO:
%    Make sure QBX knows which side of the domain the box is, do 2-sided case.
%
% Barnett 9/25/12

s = b.seg;
nv = 10; l = (1:nv)/nv;  % polygonal approx to box, nv verts on each of 4 sides
hx = h(1); hy = h(1); if numel(h)==2, hy = h(2); end
hy = 1i*hy;   % make go in imag direction
spill = 0.1*hy;   % how much to spill down below the real line to make sure
                  % all points collected even if inpolygon of s not good
b = [2*hx*l-spill 2*hx-spill+(hy+spill)*l 2*hx+hy-2*hx*l hy-(hy+spill)*l].' - hx; % vertex list in t-plane relative to center parameter value t
box.v = s.Z(t + b);

if nargout>1    % want inside-box testing... (returns [] if p=[])
 i = utils.inpolywrapper(p, box.v);  
end
