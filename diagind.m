function i = diagind(A, d)
% DIAGIND - diagonal (or sub- or super-diagonal) indices of general matrix.
%
% return diagonal indices of a general matrix, useful for changing a diagonal
% in O(N) effort, rather than O(N^2) if add a matrix to A using matlab diag()
%
% i = diagind(A) returns diagonal indices of square or rectangular matrix A.
%
% i = diagind(A, d) returns indices of d'th subdiagonal (d may be negative)
%
% barnett 2/6/08, generalized to rectangular matrices 2/23/10

if nargin<2, d = 0; end
M = size(A,1); N = size(A,2); P = min(N,M);
if N==M
  if d>=0
    i = sub2ind(size(A), 1+d:P, 1:P-d);
  else
    i = sub2ind(size(A), 1:P-abs(d), 1+abs(d):P);
  end
else                                   % rectangular case
    if d>=0
    n = min(N,M-d); i = sub2ind(size(A), d+(1:n), 1:n);
  else
    n = min(M,N-abs(d)); i = sub2ind(size(A), 1:n, abs(d)+(1:n));
  end
end