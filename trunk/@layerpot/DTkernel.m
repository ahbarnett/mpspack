function u = DTkernel(k, x, nx, y, ny)
% DTKERNEL - Kernel function for the D^T (single-layer potential n-deriv) op
%
% u = DTkernel(k, x, nx, y, ny)
% deriv of single-layer kernel k(x,y),
% without speed factor due to parametrization.
% y, ny are source location and normal vector (as C-#s), x, nx are same for
% target. All may be lists (or matrices) of same size.
% k is omega the wavenumber.
% K(s,t) = (ik/4) H_1^{(1)}(k.|x-y|). cos(angle(x-y, -n_x))
% Is adjoint of Dkernel
d = y - x; r = abs(d);
u = (1i*k/4) * besselh(1,k*r) .* real(conj(nx) .* d) ./ r;

