classdef pointset < handle

    % A simple class, which contains points plus associated normal
    % directions

    properties
        x                      % quadrature points in C (col vec)
        nx                     % outward unit normals in C (col vec)
    end

    methods
        function pts=pointset(x,nx)
            pts.x=[]; pts.nx=[];
            if nargin<2, nx=[]; end;
            if nargin<1, x=[]; end;
            pts.x=x; pts.nx=nx;
        end


    end
end
