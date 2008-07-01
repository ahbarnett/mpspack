% SEGMENT - create segment object
%
%  s = SEGMENT(M, [y0 ye]) line segment
%  s = SEGMENT(M, [yc R ti tf])     circular arc
%  s = SEGMENT(M, @(x)) analytic function on [0,1]



classdef segment <  handle & pointset
    properties
        w                      % quadrature weights, sum = seg len (row vec)
        eloc                   % [start point; end point] in C
        eang                   % [start angle; end angle] in [0,2pi)
        approxv                % vertex list for polygonal approximation
    end
    methods
        function s = segment(M, p)
            if ~isnumeric(p)
                s.x = []; s.w = []; s.nx = []; s.eloc = []; s.eang = [];
            elseif numel(p)==2
                [z s.w] = quadr.traprule(M);
                d = p(2)-p(1);
                len = abs(d);
                s.w = s.w*len/2;
                s.x = p(1) + d*(1+z)/2;
                s.nx = -1i * ones(size(s.x)) * d/len;
                s.eloc = [p(1); p(2)];
                s.eang = [1i;1i]*s.nx(1);
                s.approxv = [p(1); p(2)];   % start and endpt (one will drop)
            end
        end
        h = plot(s, pm)
    end
end

