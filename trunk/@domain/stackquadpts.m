function [x nx] = stackquadpts(s, pm)
% STACKQUADPTS - helper routine, ordered quad pts from signed connected seg list
%
% [x nx] = DOMAIN.STACKQUADPTS(segs, pm) returns list of x quadr pts in correct
%   order for a signed connected segment list, and optionally, nx normals
%   (with signs correct for the segment list).

% Copyright (C) 2008 - 2010, Alex Barnett, Timo Betcke

x = []; nx = [];
for j=1:numel(s)
  if pm(j)==1
    x = [x; s(j).x];
    if nargout>1, nx = [nx; s(j).nx]; end   % added Alex 8/13/10
  else                                      % reverse order since flipped seg
    x = [x; s(j).x(end:-1:1)];
    if nargout>1, nx = [nx; -s(j).nx(end:-1:1)]; end      % note minus sign !
  end
end
