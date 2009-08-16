% EVALBASES - evaluate all basis sets in a domain object, on a pointset
%
%  A = EVALBASES(d, p) returns matrix A whose jth column is the jth basis
%   function in domain d evaluated on the pointset p.
%
%  [A An] = EVALBASES(d, p) returns also the normal derivatives using the
%   normals associated with the pointset.
%
%  [A Ax Ay] = EVALBASES(d, p) returns A and the basis x- and y-partial
%   derivatives, ignoring the normals associated with the pointset
%
% Notes: 1) This routine is no longer used to fill problem, BVP, etc, matrices,
%   since they use ordering of basis degrees-of-freedom derived from the bas
%   objects in the problem. Basis objects in a problem do not reside inside
%   domains, so the dof ordering returned by this routine won't in general
%   match those in problem, BVP, etc. It may become obsolete.
%   2) For layer potential basis sets, and if p is also a segment
%   object, jump relations will be taken into account, corresponding to
%   evaluation on the interior (limit approaching from inside) of the domain.

% Copyright (C) 2008, 2009, Alex Barnett, Timo Betcke


function [A Ax Ay] = evalbases(d, p, opts)

A = []; Ax = []; Ay = [];   % matrices will be stacked as columns

if nargin<3, opts = []; end
for b=d.bas                    % loop over only basis set objects in domain
  bas = b{1};                  % ugly, extracts object from cell
  opts.dom = d;                % pass in which domain we're in (for jump rels)
  if nargout==1
    [bA] = bas.eval(p, opts); A = [A bA]; % stack as blocks of columns
  elseif nargout==2
    [bA bAn] = bas.eval(p, opts); A = [A bA]; Ax = [Ax bAn];
  elseif nargout==3
    [bA bAx bAy] = bas.eval(p, opts);
    A = [A bA]; Ax = [Ax bAx]; Ay = [Ay bAy]; 
  end
end
