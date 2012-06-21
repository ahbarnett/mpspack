function u = Skernel(k, x, nx, y, ny)
% SKERNEL - Kernel function for the S (single-layer potential value) operator
%
% u = Skernel(k, x, nx, y, ny)
% single-layer kernel function k(x,y),
% without speed factor due to parametrization.
% y, ny are source location and normal vector (as C-#s), x, nx are same for
% target. All may be lists (or matrices) of same size.
% k is omega the wavenumber.
% returned kernel is K(x,y) = (i/4) H_0^{(1)}(k.|x-y|)
if k==0
  u = (-1/(2*pi))*log(abs(x-y));
else
  u = (1i/4) * besselh(0,k*abs(x-y));
end
