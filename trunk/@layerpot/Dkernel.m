function u = Dkernel(k, x, nx, y, ny)
% DKERNEL - Kernel function for the D (double-layer potential value) operator
%
% u = Dkernel(k, x, nx, y, ny)
% double-layer kernel function k(x,y),
% without speed factor due to parametrization.
% y, ny are source location and normal vector (as C-#s), x, nx are same for
% target. All may be lists (or matrices) of same size.
% k is omega the wavenumber.
% K(s,t) = (ik/4) H_1^{(1)}(k.|x-y|). cos(angle(x-y, n_y))
d = x - y; r = abs(d);
if k>0
  u = (1i*k/4) * besselh(1,k*r) .* real(conj(ny) .* d) ./ r; % real(..)=dot prod
else
  u = (1/2/pi) * real(conj(ny) .* d) ./ r.^2;
end
