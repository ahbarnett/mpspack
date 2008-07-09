function sc = Jscalefactor(l,k)
sc = besselj(l, k*1.0);            % use typical max radius of 0.7 - critical
% get roundoff problems if make typ radius too small; keep >=1.0