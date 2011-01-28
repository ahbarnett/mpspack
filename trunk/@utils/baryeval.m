function u = baryeval(x, w, y, t)
%  u = baryeval(x, w, y, t) evaluates at the list of values t, the barycentric
%   interpolation formula using nodes x and weights w (found using
%   baryweights(x)), and with data y=f(x). Algorithm is O(NM)
%
% Based upon Berrut-Trefethen SIREV 2004, formula (4.2), 2nd true barycentric
%
% See also: BARYWEIGHTS, BARYPROJS

% Copyright (C) 2011 Alex Barnett

L = utils.baryprojs(x, w, t); % matrix
u = L * y(:);                 % gives a col vec
u = reshape(u, size(t));
