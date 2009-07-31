function s = polyseglist(M, p, qtype)
% POLYSEGLIST - create closed list of segment objects from CCW polygon vertices
%
%  s = POLYSEGLIST(M, p) creates segment objects with s a list of pointers to
%   them, corresponding to a closed polygon with vertices in the list p. An
%   equal number M of quadrature points are used per edge. If sense is CCW then
%   normals of segments point outwards.
%
%  s = POLYSEGLIST(M, p, qtype) uses specified quadrature type (see SEGMENT).
%
% Copyright (C) 2008, 2009, Timo Betcke, Alex Barnett

p = reshape(p, [1 numel(p)]);
nextp = circshift(p, [0 -1]);
for j=1:numel(p)
  if nargin==2
    s(j) = segment(M, [p(j) nextp(j)]);
  else
    s(j) = segment(M, [p(j) nextp(j)], qtype);
  end
end
