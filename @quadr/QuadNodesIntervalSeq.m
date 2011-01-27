function [Ax, Aw] = QuadNodesIntervalSeq(Apoints, AcorrL, AcorrR, h, order)
% Return quadrature nodes and weights on a piecewise linear contour given in Apoints
% ACorrL, ACorrR are the arrays the same size as Apoints, giving the endpoint correctiong to be used
% on the left and right side of each point (the first ACorrL and the last ACorrR values are ignored):
%        0: none
%        1: smooth function
%        2: square root singularity
%        3: log singularity
% order: is the order of the endpoint corrections used
% The nodes (Ax) and the weights (Aw) are returned

    NI = length(Apoints) - 1;
    Ax = [];
    Aw = [];

    for j=1:NI
	[Ax1, Aw1] = QuadNodesInterval(Apoints(j), Apoints(j+1), 0, h, AcorrR(j), AcorrL(j+1), order);
	Ax = [Ax; Ax1];
	Aw = [Aw; Aw1];
    end

%end
