function h = tsubplot(m,n,j)
% TSUBPLOT - basic subplot variant that places axes closer together and larger
%
%  h = TSUBPLOT(m, n, j) returns handle of axes created at j^th position in
%   m-by-n rectangular grid.

hfac = 0.95;                      % if 1, completely filled.
vfac = 0.9;                      % if 1, completely filled.
hgap = 0.5 * (1-hfac); vgap = 0.4 * (1-vfac);
h = axes('position',[(hgap+mod(j-1,n))/n 1-(-vgap+ceil(j/n))/m hfac/n vfac/m]);
