classdef regfbbasis < handle & basis
    
    % Basis of regular Fourier-Bessel functions.
    
    properties
        origin   % Origin of the Fourier-Bessel fct.
        realflag % Decide whether the basis is evaluated using real
                 % sine/cos or complex exponentials
    end
    
    methods
        function regfb = regfbbasis(origin,N,realflag,k)
            if nargin<4, k=NaN; end;
            if nargin<3, realflag=0; end; % By default use complex. exp.
            if nargin<2, N=20; end; % Default degree of FB fct.
            if nargin<1, origin=0; end; % Default origin is zero
            
            regfb.k=k;
            regfb.realflag=realflag;
            regfb.N=N;
            regfb.origin=origin;
        end
        function [A, A1, A2] = eval(regfb,pts)
        
            % Evaluates the basis at a given set of points
            np=length(pts.x); % Number of points
            R=abs(pts.x-origin);
            ang=angle(pts.x-origin);
            if realflag
                bes=besselj(0:N,k*R)
                c=cos(ang*(1:N));
                s=sine(ang*(1:N));
                A=[bes(:,1), bes(:,2:end).*c, bes(:,2:end).*s];
                if nargout>1,
                    bespp=[bes(:,2:end),besselj(N+1,k*R)];
                    besr=(repmat(0:N,np,1).*bes-k*repmat(R,N+1,1).*bespp)./R;
                    
           
            
    end
end