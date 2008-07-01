% static class of quadrature rules.

classdef quadr
    methods(Static)
        [x w] = traprule(N);
        [x w] = gauss(N);
        [x w] = clencurt(N);
    end
end
