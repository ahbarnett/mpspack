% POINTSET - create a pointset object with locations and normal vectors as
% complex numbers.
%
% A pointset is simple object containing a list of points in 2D, plus possibly
% associated normal directions.  It is used to store quadrature points on a
% segment, and also evaluation point lists. Coordinates are stored as 
% complex numbers.
%
%  p = POINTSET() creates an empty object
%  
%  p = POINTSET(x) where x is m-by-1 array, creates pointset with m points, where
%  the ith point has Cartesian coordinates (Re x(i), Im x(i)).
%  
%  p = POINTSET(x, nx) where x is above and nx has same size as x, creates
%  pointset with coordinates x (interpreted as above) and associated normals
%  nx (interpreted in the same way). The Euclidean lengths of the vectors in
%  nx are not required to be, nor changed to, unity.
%
%  p = POINTSET(s) where s is a segment or array of segments makes a pointset
%  by stacking the points and normals in the segments.
%
% See also: POINTSET/plot, SEGMENT which builds on POINTSET

% Copyright (C) 2008-2010, Alex Barnett, Timo Betcke


classdef pointset < handle

    properties
        x       % complex quadrature points (col vec)
        nx      % (optional) outward unit normals as complex numbers (col vec)
                               
    end

    methods
        function pts = pointset(x,nx)    % creator for pointset object
            pts.x=[]; pts.nx=[];
            if nargin>0
              if isa(x, 'segment')       % it's a segment (array) so get x,nx
                pts.x = vertcat(x.x); pts.nx = vertcat(x.nx);
              else
                pts.x = x;
                if ~isempty(x) && size(x,2)~=1, error('x must be m-by-1'); end
              end
            end
            if nargin>1 && ~isempty(nx)
              if size(nx)==size(x), pts.nx=nx;
              else error('nx must be same size as x'); end
            end
        end

        % external functions...
        h = plot(pts)
    end
end
