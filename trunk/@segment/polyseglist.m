function s = polyseglist(M, p, qtype, opts)
% POLYSEGLIST - create closed list of segment objects from CCW polygon vertices
%
%  s = POLYSEGLIST(M, p) creates segment objects with s a list of pointers to
%   them, corresponding to a closed polygon with vertices in the list p. An
%   equal number M of quadrature points are used per edge. If sense is CCW then
%   normals of segments point outwards.
%
%  s = POLYSEGLIST(M, p, qtype) or s = POLYSEGLIST(M, p, qtype, opts) uses
%   specified quadrature type (see SEGMENT).

% Copyright (C) 2008 - 2012, Alex Barnett, Timo Betcke

p = reshape(p, [1 numel(p)]);
nextp = circshift(p, [0 -1]);
for j=1:numel(p)
  if nargin==2
    s(j) = segment(M, [p(j) nextp(j)]);
  elseif nargin==3
    s(j) = segment(M, [p(j) nextp(j)], qtype);
  else
    s(j) = segment(M, [p(j) nextp(j)], qtype, opts);
  end
end
