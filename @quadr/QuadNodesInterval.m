function [Ax, Aw] = QuadNodesInterval(a, b, N, h, corra, corrb, order)
% Return the quadrature nodes and weights on an interval [a,b]
% with N points or step size approximately h (if N is set to 0)
% corra, corrb determine the endpoint corrections used at the two endpoints:
%        0: none
%        1: smooth function
%        2: square root singularity
%        3: log singularity
% order: is the order of the endpoint corrections used
% The nodes (Ax) and the weights (Aw) are returned

    if (corra == 0)
	NodesToSkipL = 0;
    elseif (corra == 1)
	[ExtraNodesL, ExtraWeightsL, NodesToSkipL] = quadr.QuadSmoothExtraPtNodes(order);
    elseif (corra == 2)
	[ExtraNodesL, ExtraWeightsL, NodesToSkipL] = quadr.QuadSqrtExtraPtNodes(order);
    elseif (corra == 3)
	[ExtraNodesL, ExtraWeightsL, NodesToSkipL] = quadr.QuadLogExtraPtNodes(order);
    else
	error(sprintf('Unknown endpoint corretion method corra=%d', corra));
    end

    if (corrb == 0)
	NodesToSkipR = 0;
    elseif (corrb == 1)
	[ExtraNodesR, ExtraWeightsR, NodesToSkipR] = quadr.QuadSmoothExtraPtNodes(order);
    elseif (corrb == 2)
	[ExtraNodesR, ExtraWeightsR, NodesToSkipR] = quadr.QuadSqrtExtraPtNodes(order);
    elseif (corrb == 3)
	[ExtraNodesR, ExtraWeightsR, NodesToSkipR] = quadr.QuadLogExtraPtNodes(order);
    else
	error(sprintf('Unknown endpoint corretion method corrb=%d', corrb));
    end

    if (N < 2)
	N = max(2, ceil(abs(b-a)/h));
    end
    
    % We have to have enough nodes to skip
    N = max(NodesToSkipL + NodesToSkipR, N);

    % Generate the regular nodes
    N1 = N-1; 
    h = (b-a)/N1;
    Ax = a + (0:N1)' * h;
    Aw = ones(size(Ax)) * h;
    Aw(1) = 0.5*h;
    Aw(length(Ax)) = 0.5*h;
    
    % Add the left endpoint corrections
    if (NodesToSkipL > 0)
	Ax = [a+ExtraNodesL*h; Ax(NodesToSkipL+1:length(Aw)) ];
	Aw = [ExtraWeightsL*h; Aw(NodesToSkipL+1:length(Aw)) ];
    end

    % Add the right endpoint corrections
    if (NodesToSkipR > 0)
	Ax = [ Ax(1:length(Ax)-NodesToSkipR); flipud(b-ExtraNodesR*h) ];
	Aw = [ Aw(1:length(Aw)-NodesToSkipR); flipud(ExtraWeightsR*h) ];
    end

%end
