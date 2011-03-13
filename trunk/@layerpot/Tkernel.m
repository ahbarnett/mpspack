function u = Tkernel(k, x, nx, y, ny)
% TKERNEL - Kernel function for the T (double-layer potential n-deriv) operator
%
% deriv of double-layer kernel k(x,y),
% without speed factor due to parametrization.
% y, ny are source location and normal vector (as C-#s), x, nx are same for
% target. All may be lists (or matrices) of same size.
% k is omega the wavenumber.
% K(s,t) = yukky stuff.
d = y - x; r = abs(d);
csrx = conj(nx).*d;                     % (code taken from above)
csry = conj(ny).*d;             % cos src normals
cc = real(csry).*real(csrx) ./ (r.*r);      % cos phi cos th
cdor = real(csry.*csrx) ./ (r.*r.*r);   % cos(phi-th) / r
u = (1i*k/4)*besselh(1,k*r) .* (-cdor) + (1i*k*k/4)*cc.*besselh(0,k*r);
