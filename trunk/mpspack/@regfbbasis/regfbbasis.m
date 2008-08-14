classdef regfbbasis < handle & basis

    % REGFBBASIS - create a regular Fourier-Bessel (cylindrical J-exp) basis set
    %
    %  b = REGFBBASIS(origin, N, k, opts) creates a regular FB basis
    %   object, with origin, N being max order, wavenumber k, and options:
    %   opts.real: if true, use real values (cos and sin type), otherwise exp.
    %
    % To Do:
    %
    %  - Get derivatives at zero right
    properties
        origin   % Origin of the Fourier-Bessel fct.
        realflag % Decide whether the basis is evaluated using real
                                   % sine/cos or complex exponentials
        usegsl   % Use GSL Bessel function if true
        rescale_arg  % argument (kR) to use to rescale the J bessels.
    end

    methods
        function regfb = regfbbasis(origin, N, k, opts)     % constructor
            if nargin<3, k=NaN; end;
            if nargin<4, opts = []; end
            if ~isfield(opts,'real'), opts.real = 1; end   % default
            if ~isfield(opts,'usegsl'), opts.usegsl=0; end
            if ~isfield(opts,'rescale_rad'), opts.rescale_rad=0; end
            if nargin<2, N=20; end; % Default degree of FB fct.
            if nargin<1, origin=0; end; % Default origin is zero

            regfb.k=k;
            regfb.realflag=opts.real;
            regfb.usegsl=opts.usegsl;
            regfb.N=N;
            regfb.origin=origin;
            regfb.Nf = 2*N+1;           % there are 2N+1 functions
            regfb.rescale_arg = k * opts.rescale_rad;
        end
        
        function [A, A1, A2] = eval(regfb,pts,opts) % ........... evaluator

        % Evaluates the basis at a given set of points
            resc = (regfb.rescale_arg>0);
            N=regfb.N; k=regfb.k;
            np=length(pts.x); % Number of points
            R=abs(pts.x-regfb.origin);
            ang=angle(pts.x-regfb.origin);
            [bes,err]=regfb.besselwrapper(N+1,k*R); % Use GSL function
            if nnz(err)>0,
                warning('Error in computing regular Bessel functions. Try to reduce basis size.');
            end
            c=cos(ang*(1:N));
            s=sin(ang*(1:N));
            A=[bes(:,1), bes(:,2:end-1).*c, bes(:,2:end-1).*s];
            if resc              % rescale by orders
              scfac = 1./besselj(0:N, min(0:N, regfb.rescale_arg));
              A = A .* repmat([scfac scfac(2:end)], [numel(ang) 1]);
            end
            if nargout>1, % derivs wanted
                if numel(find(R==0))>0,
                    warning('Computing x/y or normal derivatives of regular Bessel functions at origin not implemented');
                end
                
                if nnz(err)>0,
                    warning('Error in computing regular Bessel functions. Try to reduce basis size.');
                end
                besr=k/2*([-bes(:,2),bes(:,1:end-2)]-bes(:,2:end));
                bestc=-repmat(0:N,np,1).*bes(:,1:end-1).*sin(ang*(0:N));
                bests=repmat(1:N,np,1).*bes(:,2:end-1).*cos(ang*(1:N));
                Ar=[besr(:,1),besr(:,2:end).*c,besr(:,2:end).*s];
                At=[bestc,bests];
                cc=repmat(cos(ang),1,2*N+1); ss=repmat(sin(ang),1,2*N+1);
                RR=repmat(R,1,2*N+1);
                if nargout==2,
                    nx=repmat(real(pts.nx),1,2*N+1); ny=repmat(imag(pts.nx),1,2*N+1);
                    A1=Ar.*(nx.*cc+ny.*ss)+At.*(ny.*cc-nx.*ss)./RR;
                    if resc              % rescale by orders
                      A1 = A1 .* repmat([scfac scfac(2:end)], [numel(ang) 1]);
                    end
                end
                if nargout==3,
                    A1=cc.*Ar-ss.*At./RR; A2=ss.*Ar+cc.*At./RR;
                    if resc              % rescale by orders
                      A1 = A1 .* repmat([scfac scfac(2:end)], [numel(ang) 1]);
                      A2 = A2 .* repmat([scfac scfac(2:end)], [numel(ang) 1]);
                    end
                end
            end
            if ~regfb.realflag,
                sp=repmat((-1).^(N:-1:1),np,1);
                A=[sp.*(A(:,N+1:-1:2)-1i*A(:,end:-1:N+2)),A(:,1),A(:,2:N+1)+1i*A(:,N+2:end)];
                if nargout==2,
                    A1=[sp.*(A1(:,N+1:-1:2)-1i*A1(:,end:-1:N+2)),A1(:,1),A1(:,2:N+1)+1i*A1(:,N+2:end)];
                end
                if nargout==3,
                    A1=[sp.*(A1(:,N+1:-1:2)-1i*A1(:,end:-1:N+2)),A1(:,1),A1(:,2:N+1)+1i*A1(:,N+2:end)];
                    A2=[sp.*(A2(:,N+1:-1:2)-1i*A2(:,end:-1:N+2)),A2(:,1),A2(:,2:N+1)+1i*A2(:,N+2:end)];
                end
            end
        end   % ....................... end function eval
        
        function showgeom(regfb, opts) % ...................... show basis geom
          if nargin<2, opts = []; end
          plot(real(regfb.origin), imag(regfb.origin), 'ro');
          if isfield(opts, 'label')
            text(real(regfb.origin), imag(regfb.origin), opts.label);
          end
        end % func
        
    end % ... end methods
    methods (Access=private)
        function [ret,err]=besselwrapper(regfb,N,r)
            % Depending on regfb.usegsl return either Bessel's up to order
            % N computed by GSL or use Matlab implementation of recursive
            % Bessel functions

            if regfb.usegsl,
                [ret,err]=utils.gslbesselj(0,N,r);
            else
                ret=utils.recursivebessel(N,r);
                err=0; % No error checking implemented for utils.recursivebessel
            end
        end
    end
end