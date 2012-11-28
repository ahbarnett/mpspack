classdef regfbbasis < handle & basis

    % REGFBBASIS - create a regular Fourier-Bessel basis set
    %
    %  b = REGFBBASIS(origin, N) creates a regular Fourier-Bessel basis object
    %   with given origin, and order N. As with all basis types, the wavenumber
    %   k is determined by that of the affected domain, as follows:
    %
    %   k=0 gives harmonic polynomials,
    %       {1, Re z, Re z^2, ..., Re z^N, Im z, ..., Im z^N}   (real case)
    %    or {(-conj(z))^N, (-conj(z))^{N-1} ,...., -conj(z),
    %                                        1, z,..., z^N}     (complex case)
    %   k>0 gives generalized harmonic polynomials,
    %       {J_n(kr)cos(n.theta)} n=0,..,N and {J_n(kr)sin(n.theta)} n=1,..,N
    %    or {J_n(kr)exp(in.theta)} n=-N,..,N                    (complex case)
    %
    %  b = REGFBBASIS(origin, N, opts) is as above, with options:
    %   opts.real: if true (default), use real values (cos, sin type), otherwise
    %              use complex exponentials with orders -N through N.
    %   opts.rescale_rad: if positive, rescales the basis coefficients i.e.
    %              columns of evaluation matrix, such that the value at radius
    %              rescale_rad is O(1). (default is 0, giving no rescaling)
    %   opts.besselcode: math library to use for J Bessel evaluation (k>0)
    %              = 'r' use downwards-recurrence in Matlab, fast, but
    %                        relative accuracy not guaranteed.
    %              = 'm' use Matlab's built-in besselj, is slower (default).
    %              = 'g' use GNU Scientific Library via MEX interface, fast.
    %
    % Note: the alternating signs in the negative orders for k=0, real=0 case
    % is to match the signs of the k>0 Bessels reflection principles.
    %
    % See also: DOMAIN/ADDREGFBBASIS
    
    % Copyright (C) 2008 - 2012, Alex Barnett, Timo Betcke
    
    
    properties
        origin             % Origin of the Fourier-Bessel fct.
        real               % true for real valued output, otherwise complex
        rescale_rad        % radius to use to rescale the polys or J bessels.
                           %  (changed from rescale_arg by barnett 7/14/09)
        besselcode         % character controlling special function library used
        fastang            % if true, use recurrent trig funcs (not implemented)
        ViS         % experimental, preconditions coeffs, ie postmultiply A, etc
    end

    methods
        function regfb = regfbbasis(origin, N, opts)     % constructor
            if nargin<3, opts = []; end
            if nargin<2 | isempty(N), N=20; end;          % default degree
            if nargin<1 | isempty(origin), origin=0; end; % default origin
            if ~isfield(opts,'real'), opts.real = 1; end  % default opts...
            if ~isfield(opts,'rescale_rad'), opts.rescale_rad=0; end
            if ~isfield(opts,'besselcode'), opts.besselcode='m'; end
            if ~isfield(opts,'fastang'), opts.fastang=0; end
            if ~isfield(opts,'nmultiplier'), opts.nmultiplier=1; end

            regfb.nmultiplier=opts.nmultiplier;
            regfb.N=N*regfb.nmultiplier; 
            regfb.origin=origin;
            regfb.real=opts.real; regfb.besselcode=opts.besselcode;
            regfb.rescale_rad = opts.rescale_rad;
            regfb.fastang = opts.fastang;
            
            
            if isempty(strfind('rmg', regfb.besselcode))  % usage check
              error('opts.besselcode must be one of ''r'',''m'',''g''!');
            end
        end
        
        function Nf = Nf(regfb) % ...................... returns # dofs
          Nf = 2*regfb.N + 1;
        end
        
        function [A, A1, A2] = eval(regfb,pts,opts) % ........... evaluator
        % EVAL - evaluates regular Fourier-Bessel basis at given set of points
        %
        % A = EVAL(p) where p is a pointset object containing M points, returns
        %   a M-by-Nf matrix whose ij'th entry is Phi_j(z_i), where Phi_j is
        %   the jth basis function, and z_i the ith point. Nf is the number
        %   of degrees of freedom in the basis set object.
        %        
        % [A An] = EVAL(p) also returns matrix An whose ij'th entry is
        %   d/dn_i Phi_j(z_i), the derivative in the ith normal direction in
        %   the pointset.
        %
        % [A Ax Ay] = EVAL(p) returns A as above, and matrices of basis function
        %    derivatives in the x- and y-directions. That is, Ax has ij'th
        %    entry d/dx Phi_j(z_i) while Ay has ij'th entry d/dy Phi_j(z_i)
        %
        % Also see: POINTSET, REGFBBASIS
            resc = (regfb.rescale_rad>0); % whether to rescale or not
            N = regfb.N;
            k = regfb.k;                  % NB this is now a method not property
            if isnan(k), error('regfbbasis: no domain wavenumber defined!'); end
            np = length(pts.x);               % Number of eval points
            R=abs(pts.x-regfb.origin);        % col vec of distances
            ang=angle(pts.x-regfb.origin);    % col vec of angles
            if k==0                % Laplace (harmonic polynomials 0,1,..,N+1)
              bes = repmat(R, [1 N+2]).^repmat(0:N+1, [numel(R) 1]); % eval R^n
            else                   % Helmholtz J_n for n=0...N+1 (for derivs N)
              [bes, err]=regfb.besselwrapper(N+1, k*R); % eval w/ a bessel code
              if nnz(err)>0
                warning('Error in computing regular Bessel functions. Try to reduce basis size.');
              end
            end
            if regfb.fastang      % intelligent exps using complex rotations
              %angs = utils.fast_angle_exp(0:N+1, pts.x-regfb.origin); TODO
              error('fastang not yet implemented');
            else                  % conventional eval
              c=cos(ang*(1:N));   % NB * here is outer product! (also below)
              s=sin(ang*(1:N));
            end
            % stack cos funcs 0...N then sin funcs 1...N :
            A = [bes(:,1), bes(:,2:end-1).*c, bes(:,2:end-1).*s];
            if resc                                    % rescale by orders
              scfac = 1./regfb.Jrescalefactors(0:N);   % uses rescale_arg
              A = A .* repmat([scfac scfac(2:end)], [np 1]);
            end
            
            if nargout>1         % --------- derivs wanted ---------------
              i = find(R==0); R(i) = 1;     % zero-R indices, use dummy R vals
              if k==0        % evaluate radial deriv n.R^(n-1) for n=0...N
                besr = [0*R repmat(1:N, [np 1]).*bes(:,1:end-2)];
              else
                besr = k/2*([-bes(:,2),bes(:,1:end-2)]-bes(:,2:end)); % dJ/dr
              end
              % Ar is derivative of radial func, At tangential derivative
              Ar = [besr(:,1), besr(:,2:end).*c, besr(:,2:end).*s];
              clear besr;
              At = [-repmat(0:N,np,1).*bes(:,1:end-1) .* [0*ang, s], ...
                    repmat(1:N,np,1).*bes(:,2:end-1) .* c];  % reuse s & c
              clear c s bes;
              cc=repmat(cos(ang),1,2*N+1); ss=repmat(sin(ang),1,2*N+1);
              invRR=repmat(1./R,1,2*N+1);  % saves division later which is slow
              
              if nargout==2      % normal directional deriv ------
                nx=repmat(real(pts.nx),1,2*N+1);
                ny=repmat(imag(pts.nx),1,2*N+1);
                A1 = Ar.*(nx.*cc+ny.*ss) + At.*(ny.*cc-nx.*ss).*invRR;
                if ~isempty(i)  % overwrite the zero-R derivs (why k/2?)
                  A1(i,:) = 0;
                  if k==0,
                    A1(i,2) = nx(i,1); % n=1 cos term
                    A1(i,N+2) = ny(i,1);            % n=1 sin term
                  else
                    A1(i,2) = k/2*nx(i,1); % n=1 cos term
                    A1(i,N+2) = k/2*ny(i,1);            % n=1 sin term
                  end
                end
                if resc              % rescale by orders
                  A1 = A1 .* repmat([scfac scfac(2:end)], [np 1]);
                end
                
              elseif nargout==3     % x- and y-derivs --------
                A1 = cc.*Ar - ss.*At.*invRR; A2 = ss.*Ar + cc.*At.*invRR;
                if ~isempty(i)  % overwrite the zero-R derivs (why k/2?)
                  A1(i,:) = 0; A2(i,:) = 0;
                  if k==0, A1(i,2) = 1; A2(i,N+2) = 1;
                  else, A1(i,2) = k/2; A2(i,N+2) = k/2; end
                end
                if resc              % rescale by orders
                  A1 = A1 .* repmat([scfac scfac(2:end)], [np 1]);
                  A2 = A2 .* repmat([scfac scfac(2:end)], [np 1]);
                end
              end
            end

            if ~regfb.real   % construct complex combos from sin,cos
              %  this works even for the zero-radius points
              sp=repmat((-1).^(N:-1:1),np,1);
              A=[sp.*(A(:,N+1:-1:2)-1i*A(:,end:-1:N+2)),A(:,1),...
                 A(:,2:N+1)+1i*A(:,N+2:end)];
              if nargout==2
                A1=[sp.*(A1(:,N+1:-1:2)-1i*A1(:,end:-1:N+2)),A1(:,1),...
                    A1(:,2:N+1)+1i*A1(:,N+2:end)];
              elseif nargout==3
                A1=[sp.*(A1(:,N+1:-1:2)-1i*A1(:,end:-1:N+2)),A1(:,1),...
                    A1(:,2:N+1)+1i*A1(:,N+2:end)];
                A2=[sp.*(A2(:,N+1:-1:2)-1i*A2(:,end:-1:N+2)),A2(:,1),...
                    A2(:,2:N+1)+1i*A2(:,N+2:end)];
              end
            end

            if ~isempty(regfb.ViS), A = A*regfb.ViS; % precond, experimental...
              if nargout>1, A1 = A1*regfb.ViS; end 
              if nargout>2, A2 = A2*regfb.ViS; end
            end                                      % ...end experimental!

        end   % ....................... end function eval
        
        function showgeom(regfb, opts) % ...................... show basis geom
        % SHOWGEOM - plot regular FB basis set geometry info (just origin now)
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
            sc = regfb.rescale_rad.^n;
          else                      % Bessel case
            sc = abs(besselj(n, min(n, k*regfb.rescale_rad))); % n can be list.
            % The min here stops the J from getting to osc region (it stays 
            % at turning point). Note Matlab's bessel is used for relative
            % accuracy when the value is very small.
          end
        end
        
    end % ... end methods
    
    methods (Access=private)   % user can't see ----------------------------
        function [ret,err]=besselwrapper(regfb,N,r)
        % Interface to J-Bessel special function libraries (see opts in
        % regfbbasis). Returns matrix of J_0 through J_N evaluated at the
        % argument list in r. Matrix size is M-by-N+1, where M=numel(r).
          switch regfb.besselcode
           case 'r'
            ret=utils.recurrencebesselJ(N,r);    % Alex's recurrence code
            err=0; % No error checking implemented for utils.recurrencebesselJ
           case 'g'
            [ret,err]=utils.gslbesselj(0,N,r); % Timo's MEX interface to GSL
           case 'm'
            ret = besselj(repmat(0:N, [numel(r) 1]), repmat(r(:), [1 N+1]));
            %ret = bsxfun(@besselj, 0:N, r(:)); % slow, calls besschk loads
            err = 0;   % >2012a releases have no err output.
           otherwise
            error('regfbbasis: invalid besselcode property!');
          end
        end % end func
    end % private methods
end