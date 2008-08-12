% EVALBASES - evaluate all basis sets in a domain object, on a pointset
%
%  A = EVALBASES(d, p) returns matrix A whose jth column is the jth basis
%   function in domain d evaluated on the pointset p.
%
%  [A An] = EVALBASES(d, p) returns also the normal derivatives using the
%   normals associated with the pointset
%
%  [A Ax Ay] = EVALBASES(d, p) returns A and the basis x- and y-partial
%   derivatives, ignoring the normals associated with the pointset
%
%  Note that for layer potential basis sets, and if p is also a segment
%  object, jump relations will be taken into account, corresponding to
%  evaluation on the interior (limit approaching from inside) of the domain.

function [A Ax Ay] = evalbases(d, p)

A = []; Ax = []; Ay = [];   % matrices will be stacked as columns

opts = [];
for b=d.bas                    % loop over basis set objects in domain
  bas = b{1};                  % ugly, extracts object from cell
  segind = find(d.seg==p);     % if p is part of domain bdry, tell which side
  if numel(segind)==1
    opts.layerpotside = -d.pm(segind);  % tell eval to take limit on inside
  elseif numel(segind)>1
    fprintf('warning: weirdness, p segment matches >1 domain segments!\n');
  end
  if nargout==1
    [bA] = bas.eval(p, opts); A = [A bA]; % stack as blocks of columns
  elseif nargout==2
    [bA bAn] = bas.eval(p, opts); A = [A bA]; Ax = [Ax bAn];
  elseif nargout==3
    [bA bAx bAy] = bas.eval(p, opts);
    A = [A bA]; Ax = [Ax bAx]; Ay = [Ay bAy]; 
  end
end
