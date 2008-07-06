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
            if nargin<3, realflag=1; end; % By default use real arithmetic
            if nargin<2, N=20; end; % Default degree of FB fct.
            if nargin<1, origin=0; end; % Default origin is zero

            regfb.k=k;
            regfb.realflag=realflag;
            regfb.N=N;
            regfb.origin=origin;
        end
        function [A, An, Ax, Ay] = eval(regfb,pts)

            % Evaluates the basis at a given set of points
            N=regfb.N; k=regfb.k;
            np=length(pts.x); % Number of points
            R=abs(pts.x-regfb.origin);
            ang=angle(pts.x-regfb.origin);
            bes=besselj(0:N,k*R);
            c=cos(ang*(1:N));
            s=sin(ang*(1:N));
            A=[bes(:,1), bes(:,2:end).*c, bes(:,2:end).*s];
            if nargout>1,
                bespp=[bes(:,2:end),besselj(N+1,k*R)];
                besr=(repmat(0:N,np,1).*bes-k*repmat(R,1,N+1).*bespp)./repmat(R,1,N+1);
                bestc=-repmat(0:N,np,1).*bes.*sin(ang*(0:N));
                bests=repmat(1:N,np,1).*bes(:,2:end).*cos(ang*(1:N));
                Ar=[besr(:,1),besr(:,2:end).*c,besr(:,2:end).*s];
                At=[bestc,bests];
                Ax=repmat(cos(ang),1,2*N+1).*Ar-repmat(sin(ang),1,2*N+1).*At./repmat(R,1,2*N+1);
                Ay=repmat(sin(ang),1,2*N+1).*Ar+repmat(cos(ang),1,2*N+1).*At./repmat(R,1,2*N+1);
                nx1=real(pts.nx); nx2=imag(pts.nx);
                An=repmat(nx1,1,2*N+1).*Ax+repmat(nx2,1,2*N+1).*Ay;
            end
            if ~regfb.realflag,
                A=[A(:,2:N+1)-1i*A(:,N+2:end),A(:,1),A(:,2:N+1)+1i*A(:,N+2:end)];
                if nargout>1,
                    An=[An(:,2:N+1)-1i*An(:,N+2:end),An(:,1),An(:,2:N+1)+1i*An(:,N+2:end)];
                    Ax=[Ax(:,2:N+1)-1i*Ax(:,N+2:end),Ax(:,1),Ax(:,2:N+1)+1i*Ax(:,N+2:end)];
                    Ay=[Ay(:,2:N+1)-1i*Ay(:,N+2:end),Ay(:,1),Ay(:,2:N+1)+1i*Ay(:,N+2:end)];
                end
            end
        end
    end
end