function w = baryweights(x)
%  w = baryweights(x) computes barycentric Lagrange weights for interpolation
%   x is an input vector of interpolation nodes, w is output of weights, same
%   size as x. The algorithm is O(N^2). No attempt to prevent over/underflow.
%
%  Based on Berrut-Trefethen SIREV 2004 paper
%
% See also: BARYEVAL

%  Copyright (C) 2011 Alex Barnett
x = x(:);   % make a column vector
N = numel(x);
X = repmat(x, [1 N]);
w = 1./prod(X - X.' + eye(N), 1);
w = reshape(w, size(x));
