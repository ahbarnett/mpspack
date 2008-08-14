classdef utils

  methods(Static)
    % things which evaluate special functions...
    [res,err]=gslbesselj(nmin,nmax,x);
    [res,err]=gslbesseljnu(v,x);
    J=recursivebessel(M,x);
    A = fundsol(r, k)
    [B radderivs] = fundsol_deriv(r, cosphi, k, radderivs)
    % other helper routines...
    b = copy(a)
    monochrome(h, c)
    u = unique(c)
    i = isin(b, c)
  end

end
