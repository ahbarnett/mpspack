classdef utils

  methods(Static)
    [res,err]=gslbesselj(nmin,nmax,x);
    [res,err]=gslbesseljnu(v,x);
  end

end
