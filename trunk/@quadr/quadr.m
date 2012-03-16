% static class of quadrature rules.

% Copyright (C) 2008 - 2011, Alex Barnett, Timo Betcke

classdef quadr
    methods(Static)
        [x w] = peritrap(N)
        [x w] = traprule(N)
        [x w] = gauss(N)
        [x w] = clencurt(N)
        [x,w,cs,ier]=kapurtrap(n,m)
        Rjn = kress_Rjn(n)
        D = perispecdiffrow(N)
        [g] = interptrig(f, N);
        
        % the following are Alpert quadrature endpoint correction rules
        % (from Andras Pataki):  test_Alpert_Pataki.m to test
        [ExtraNodes, ExtraWeights, NodesToSkip] = QuadLogExtraPtNodes(order)
        [Ax, Aw] = QuadNodesInterval(a, b, N, h, corra, corrb, order)
        [Ax, Aw] = QuadNodesIntervalSeq(Apoints, AcorrL, AcorrR, h, order)
        [ExtraNodes, ExtraWeights, NodesToSkip] = QuadSmoothExtraPtNodes(order)
        [ExtraNodes, ExtraWeights, NodesToSkip] = QuadSqrtExtraPtNodes(order)
        A = alpertizeselfmatrix(A, k, s, kerfun, opts)
        
        % roll-off functions
        f = smoothedstep(t, opts)
    end
end
