classdef utils
% UTILS - class of utility routines for mpspack.

% Copyright (C) 2008 - 2011, Timo Betcke, Alex Barnett

  methods(Static)
    
    % math libraries which evaluate special functions...
    [res,err]=gslbesselj(nmin,nmax,x);
    [res,err]=gslbesseljnu(v,x);
    J=recurrencebesselJ(M,x);
    [F0 F1 F2] = fundsol(r, k, orders, fast)
    [B radderivs] = fundsol_deriv(r, cosphi, k, radderivs)
    [H0 H1] = greengardrokhlinhank103(z);
    [H0 H1] = greengardrokhlinhank106(z);
    
    % Greengard-Gimbutas HFMM2D library...
    [U]=hfmm2dparttarg(iprec,zk,nsource,source,ifcharge,charge,ifdipole,...
                       dipstr,dipvec,ifpot,iffld,ifhess,ntarget,target,...
                       ifpottarg,iffldtarg,ifhesstarg);
    [pot,fld,hess,ier]=hfmm2dpart(iprec,zk,nsource,source,ifcharge,charge,...
                                  ifdipole,dipstr,dipvec);

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
    [f0 err N] =  extrap(f, hmax, opts)
    c = goodcaxis(u)
    h = arrow(x, y, varargin)
    
    % rootfinding and linear algebra helpers...
    [x e y u ier] = intervalrootsboyd(f, int, o)
    [r e] = trigpolyzeros(F, opts)
    [u s v info] = minsingvalvecs(A, opts)

    % interpolation...
    w = baryweights(x)
    L = baryprojs(x, w, t)
    u = baryeval(x, w, y, t)
  end

end
