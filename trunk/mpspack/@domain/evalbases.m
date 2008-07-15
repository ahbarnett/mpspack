% EVALBASES - evaluate all basis sets in a domain object, on a pointset
%
%   Issues:
%   * it's ugly that each b is a 1x1 cell which then needs a {1} reference

function [A A1 A2] = evalbases(d, p)

A = []; An = []; Ax = []; Ay = [];   % matrices will be stacked as columns

for b=d.bas                          % loop over basis set objects in domain
  if nargout==1
    [bA] = b{1}.eval(p); A = [A bA]; % stack as blocks of columns
  elseif nargout==2
    [bA bA1] = b{1}.eval(p); A = [A bA]; A1 = [A1 bA1];
  else
    [bA bA1 bA2] = b{1}.eval(p);
    A = [A bA]; A1 = [A1 bA1]; A2 = [A2 bA2]; 
  end
end
