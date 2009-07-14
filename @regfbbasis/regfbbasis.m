classdef regfbbasis < handle & basis

    % REGFBBASIS - create a regular Fourier-Bessel (cylindrical J-exp) basis set
    %
    %  b = REGFBBASIS(origin, N, opts) creates a regular FB basis
    %   object, with origin, N being max order, and options:
    %   opts.real: if true, use real values (cos and sin type), otherwise exp.
    %
    % Issues:
    %  - Get derivatives at zero right
    %  * Code radial derivatives for k=0
    properties
        origin             % Origin of the Fourier-Bessel fct.
        realflag           % Decides whether the basis is evaluated using real
                                   % sine/cos or complex exponentials
        usegsl             % Use GSL Bessel function if true
        rescale_rad        % radius to use to rescale the polys or J bessels.
                           %  (changed from rescale_arg by barnett 7/14/09)
    end

    methods
        function regfb = regfbbasis(origin, N, opts)     % constructor
            if nargin<3, opts = []; end
            if nargin<2 | isempty(N), N=20; end; % default degree of FB fct.
            if nargin<1 | isempty(origin), origin=0; end; % default origin is 0
            if ~isfield(opts,'real'), opts.real = 1; end   % default opts...
            if ~isfield(opts,'usegsl'), opts.usegsl=0; end
            if ~isfield(opts,'rescale_rad'), opts.rescale_rad=0; end

            regfb.realflag=opts.real;
            regfb.usegsl=opts.usegsl;
            regfb.N=N;
            regfb.origin=origin;
            regfb.Nf = 2*N+1;           % there are 2N+1 functions
            regfb.rescale_rad = opts.rescale_rad;
        end
        
        function [A, A1, A2] = eval(regfb,pts,opts) % ........... evaluator

        % Evaluates the reg FB basis at a given set of points. TODO write
        
            resc = (regfb.rescale_rad>0); % whether to rescale or not
            N=regfb.N;
            k = regfb.k;                  % NB this is now a method not property
            np=length(pts.x); % Number of points
            R=abs(pts.x-regfb.origin);
            ang=angle(pts.x-regfb.origin);
            if k==0                % Laplace (harmonic polynomials 0,1,..,N+1)
              bes = repmat(R, [1 N+2]).^repmat(0:N+1, [numel(R) 1]); % eval R^n
            else                   % Helmholtz
              [bes, err]=regfb.besselwrapper(N+1, k*R); % eval w/ a bessel code
              if nnz(err)>0,
                warning('Error in computing regular Bessel functions. Try to reduce basis size.');
              end
            end
            c=cos(ang*(1:N));   % NB * here is outer product! (also below)
            s=sin(ang*(1:N));
            A=[bes(:,1), bes(:,2:end-1).*c, bes(:,2:end-1).*s];
            if resc              % rescale by orders
              scfac = 1./regfb.Jrescalefactors(0:N);  % uses rescale_arg
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
        
        function sc = Jrescalefactors(regfb, n) % ... used to effect rescale_rad
        % JRESCALEFACTORS - given list of orders, return FB J-rescaling factors
        %
        %  Also now works for k=0, ie harmonic polynomials.
          k = regfb.k;
          if k==0                   % polynomial case
            sc = rescale_rad.^n;
          else                      % Bessel case
            sc = abs(besselj(n, min(n, k*regfb.rescale_rad))); % n can be list
            % the min here stops the J from getting to osc region (it stays 
            % at turning point)
          end
        end
        
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