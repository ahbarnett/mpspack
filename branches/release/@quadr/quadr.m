% static class of quadrature rules.

classdef quadr
    methods(Static)
        [x w] = peritrap(N)
        [x w] = traprule(N)
        [x w] = gauss(N)
        [x w] = clencurt(N)
        [x,w,cs,ier]=kapurtrap(n,m)
        Rjn = kress_Rjn(n)
        D = perispecdiffrow(N)
    end
end
