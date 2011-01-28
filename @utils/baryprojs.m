function L = baryprojs(x, w, t)
% L = baryprojs(x, w, t) returns a M-by-N matrix whose rows are the vectors
%  against which a data vector y could be inner producted to give its
%  interpolant at each of M target locations t. The N interpolation pts are in x
%  and the N barycentric weights in w. Vectorized for multiple t.
%  In other words, j^th col of L is the j^th Lagrange basis function eval at t.
%  O(NM)
%
% Based upon Berrut-Trefethen SIREV 2004, formula (4.2), 2nd true barycentric
%
% See also: BARYWEIGHTS, BARYEVAL

% Copyright (C) 2011 Alex Barnett

M = numel(t); N = numel(x);
diffs = repmat(t(:), [1 N]) - repmat(x(:).', [M 1]);  % M-by-N matrix
wpoles = repmat(w(:).', [M 1]) ./ diffs;              % "
denom = sum(wpoles, 2);   % col vector for each t value
L = repmat(1./denom, [1 N]) .* wpoles;                % "
L(find(diffs==0)) = 1;   % fix cases where infinity in numerator and denom
