classdef utils
% UTILS - class of utility routines for mpspack.

% Copyright (C) 2008, 2009, Timo Betcke, Alex Barnett

  methods(Static)
    
    % math libraries which evaluate special functions...
    [res,err]=gslbesselj(nmin,nmax,x);
    [res,err]=gslbesseljnu(v,x);
    J=recurrencebesselJ(M,x);
    [F0 F1 F2] = fundsol(r, k, orders, fast)
    [B radderivs] = fundsol_deriv(r, cosphi, k, radderivs)
    [H0 H1] = greengardrokhlinhank103(z);
    [H0 H1] = greengardrokhlinhank106(z);
    
    % inpolygon replacements...
    [cn,on] = inpoly(p,node,edge,TOL)
    i = inpolyc(p,v)
    i = inpolywrapper(p, v)
    
    % other helper routines...
    b = copy(a)
    monochrome(h, c)
    u = unique(c)
    i = isin(b, c)
    s = merge(s1,s2)
  end

end
