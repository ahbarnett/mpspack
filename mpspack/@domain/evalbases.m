% EVALBASES - evaluate all basis sets in a domain object, on a pointset
%
%   Issues:
%   * it's ugly that each b is a 1x1 cell which then needs a {1} reference

function [A An Ax Ay] = evalbases(d, p)

A = []; An = []; Ax = []; Ay = [];   % matrices will be stacked as columns

for b=d.bas                          % loop over basis set objects in domain
  if nargout==1
    [bA] = b{1}.eval(p); A = [A bA]; % stack as blocks of columns
  elseif nargout==2
    [bA bAn] = b{1}.eval(p); A = [A bA]; An = [An bAn];
  else
    [bA bAn bAx bAy] = b{1}.eval(p);
    A = [A bA]; An = [An bAn]; Ax = [Ax bAx]; Ay = [Ay bAy]; 
  end
end
