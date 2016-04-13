function C = Cennenbach(p)
% CENNENBACH - compute const (in small u_n norm limit) for Thm 7, Ennenbach '95
%
% C = Cennenbach(p) where p is evp object containing domain with single segment
% bdry. Star-shaped case only for now.
%
% Barnett 12/6/15

s = p.segs; d = p.doms;
xn = real(conj(s.x).*s.nx);
R = max(abs(s.x)); r = min(abs(s.x)); h = min(xn);
xi2LB = (r*h/R^3)/(R^4/r/h + R^2); % lower bnd on xi_2, Stekloff, Ennen (3.5)
xnnrmsq = sum(s.w'.*xn.^2);      % terms needed in beta
o = p.gridboundingbox;   % build grid for interior norms of geom funcs
n = floor((o.bb(2)-o.bb(1))/o.dx); gx = o.bb(1) + o.dx*(0:n); % grids
n = floor((o.bb(4)-o.bb(3))/o.dx); gy = o.bb(3) + o.dx*(0:n);      
[xx yy] = meshgrid(gx, gy); zz = xx(:) + 1i*yy(:);  % keep zz rect array
di = d.inside(zz);
r2int = o.dx^2*sum(abs(zz).^2.*di);
r4int = o.dx^2*sum(abs(zz).^4.*di);
beta = (sqrt(d.perim)/2/d.area)*( 1/sqrt(xi2LB)*sqrt(xnnrmsq - 4*d.area^2/d.perim) + .5*sqrt(r4int - r2int/d.area) );   % beta defined in Ennen Lemma 6
C = sqrt(1/xi2LB + beta^2);  % around 7.4 for rfn
