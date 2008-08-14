function h = tsubplot(m,n,j)
% TSUBPLOT - basic subplot variant that places axes closer together and larger
%
%  h = TSUBPLOT(m, n, j) returns handle of axes created at j^th position in
%   m-by-n rectangular grid.

fac = 0.95;                      % if 1, completely filled.
gap = (1-fac)/2;
h = axes('position',[(gap+mod(j-1,n))/n 1-(gap+ceil(j/n))/m fac/n fac/m]);
