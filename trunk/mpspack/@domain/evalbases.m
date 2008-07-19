% EVALBASES - evaluate all basis sets in a domain object, on a pointset
%
%   Issues:
%   * it's ugly that each b is a 1x1 cell which then needs a {1} reference

function [A Ax Ay] = evalbases(d, p)

A = []; Ax = []; Ay = [];   % matrices will be stacked as columns

for b=d.bas                          % loop over basis set objects in domain
  bas = b{1};                        % ugly
  if nargout==1
    [bA] = bas.eval(p); A = [A bA]; % stack as blocks of columns
  elseif nargout==2
    [bA bAn] = bas.eval(p); A = [A bA]; Ax = [Ax bAn];
  elseif nargout==3
    [bA bAx bAy] = bas.eval(p);
    A = [A bA]; Ax = [Ax bAx]; Ay = [Ay bAy]; 
  end
end
