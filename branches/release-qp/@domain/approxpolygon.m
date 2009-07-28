function v = approxpolygon(s, pm)
% APPROXPOLYGON - extracts col vec of approx vertices from closed seg array
%
% APPROXPOLYGON - stack together approximating polygon vertices of segment list
%
%  v = APPROXPOLYGON(seg, pm) returns a column vector of vertices (as
%   C-numbers) for the approximating polygon of a segment list seg and sign
%   list pm. If the segment list is not closed, v should not be interpreted
%   as a closed polygon.
%
% Copyright (C) 2008, 2009, Timo Betcke, Alex Barnett


v = [];
for j=1:length(s)
  if pm(j)==1
    v = [v; s(j).approxv(1:end-1)];          % drop the last point
  else
     v = [v; s(j).approxv(end:-1:2)];
  end
end
