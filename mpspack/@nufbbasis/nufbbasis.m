classdef nufbbasis < handle & basis

    % Basis of irregular Fourier-Bessel functions.

    properties
        origin   % Origin of the Fourier-Bessel fct.
        nu       % Real parameter that determines the fractional order of
                 % the Bessel fct.
        branch   % Complex angle that determines the branch cut
        offset   % Complex angle that determines the zero line of the Bessel
                 % fct.
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
                A=[A(:,1:N)-1i*A(:,N+1:end),A(:,0),A(:,1:N)+1i*A(:,N+1:end)];
                if nargout>1,
                    An=[An(:,1:N)-1i*An(:,N+1:end),An(:,0),An(:,1:N)+1i*An(:,N+1:end)];
                    Ax=[Ax(:,1:N)-1i*Ax(:,N+1:end),Ax(:,0),Ax(:,1:N)+1i*Ax(:,N+1:end)];
                    Ay=[Ay(:,1:N)-1i*Ay(:,N+1:end),Ay(:,0),Ay(:,1:N)+1i*Ay(:,N+1:end)];
                end
            end
        end
    end
end