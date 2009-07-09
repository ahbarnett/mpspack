% POINTSET - create a pointset object with locations and normal vectors in C
%
% A simple class, which contains points plus associated normal
% directions.
%
% p = POINTSET() creates empty object
% p = POINTSET(x) where x is m-by-1 list of values in C creates pointset with
%   coords x
% p = POINTSET(x, nx) where x is above and nx has same size as x, creates
%   pointset with coords x and associated normals nx. The length of normals
%   is not enforced to be 1 (this allows hacking to take various derivatives)
% 
% Notes: extra size checks added, then removed (!), Barnett 2/25/09, 5/1/09
%
% See also: POINTSET/plot, SEGMENT which builds on POINTSET

classdef pointset < handle

    properties
        x                      % quadrature points in C (col vec)
        nx                     % (optional) outward unit normals in C (col vec)
    end

    methods
        function pts = pointset(x,nx)    % creator for pointset object
            pts.x=[]; pts.nx=[];
            if nargin>0
              pts.x = x;
              % if size(x,2)~=1, error('x must be m-by-1'); end
              % error removed since other sizes useful in eg gridsolution
            end
            if nargin>1 & ~isempty(nx)
              if size(nx)==size(x), pts.nx=nx;
              else error('nx must be same size as x'); end
            end
        end

        % external functions...
        h = plot(pts)
    end
end
