classdef mfsbasis < handle & basis
% MFSBASIS - create a fundamental solutions (MFS) basis set
%
%  b = MFSBASIS(Z, N, opts) creates an MFS basis using charge points
%  at Z(t), where t are equally distributed in [0,1]. Z is a handle
%  to a given [0,1]-periodic analytic function.
%
%  b = MFSBASIS({Z, Zp}, N, opts) is as above. But Zp additionally specifies
%  the derivative of Z. This is necessary if double layer potentials are
%  used.
%
%  b = MFSBASIS(p, opts) describes MFS charge points in terms of a point
%  set object p. If double layer potentials are used then p also needs to
%  contain normal directions.
%
%  opts is an optional structure containing optional fields:
%  opts.eta   - (inf) If eta=inf use single layer potential MFS basis. For
%               eta neq inf use linear combination of double and single layer
%               potential MFS basis
%  opts.fast  - (0) Hankel evaluation method (0=matlab; 1,2=faster)
%  opts.real  - (false) If true, use real part of Hankel functions only 
properties
    Z                     % function handle for charge curve
    Zp                    % function handle for derivative
    t                     % row vec of parameter values in [0,1]
    y                     % row vec of charge points (note not col vec)
    ny                    % row vec of normal directions at charge points
    realflag              % if true, use real part only, ie Y_0
    fast                  % Hankel evaluation method (0=matlab; 1=faster)
    eta                   % If eta=inf (default) use single layer potential 
                          % MFS basis. For eta\neq inf
                          % use linear combination of double and single
                          % layer potential MFS basis
  end
  
  methods
      function b = mfsbasis(varargin) % Constructor
          if nargin==0, return; end
          if nargin>3, error('Two many arguments'); end
          defaultopts=struct('eta',inf,'fast',2,'real',0,'nmultiplier',1);
          % Now evaluate the options array
          if isstruct(varargin{end}),
              opts=utils.merge(varargin{end},defaultopts);
          else
              opts=defaultopts;
          end
          if opts.fast==2 & exist('@utils/greengardrokhlinhank106.mexglx')~=3
            opts.fast = 1;    % downgrade the speed 2->1 if 106 not available
          end
          if opts.fast==1 & exist('@utils/greengardrokhlinhank103.mexglx')~=3
            opts.fast = 0;    % downgrade the speed 1->0 if 103 not available
          end
          b.fast=opts.fast;
          b.realflag=opts.real;
          b.eta=opts.eta;
          b.nmultiplier=opts.nmultiplier;

          % Evaluate the input arguments
          
          if isa(varargin{1},'pointset'),
              % First argument pointset
              if nargin>1 && ~isstruct(varargin{2}), 
                  error('Wrong arguments'); 
              end
              b.y=(varargin{1}.x).';
              b.ny=(varargin{1}.nx).';
              if ~isempty(b.ny), b.ny=b.ny./abs(b.ny); end
              if isempty(b.ny) && ~isinf(b.eta),
                  error('Normal directions not supplied');
              end
          elseif iscell(varargin{1}),
              % First argument is cell array
              if length(varargin{1})~=2, error('Wrong argument'); end
              b.Z=varargin{1}{1};
              b.Zp=varargin{1}{2}; 
              b.updateN(varargin{2});   % set up charge pts and normals
          elseif isa(varargin{1},'function_handle'),
              % First argument is function handle
              b.Z=varargin{1};
              b.updateN(varargin{2});   % set up charge pts
              if ~isinf(b.eta), 
                  error('Derivative of Z needed for eta<inf');
              end
          else
              error('Wrong arguments');
          end
      end
      
      function [A B C] = eval(b, pts, opts)
        % EVAL - evaluate fundamental solutions (MFS) basis at set of points
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
        % Also see: POINTSET, MFSBASIS      
      o=struct('fast',b.fast);
     
      if nargin==3,                   % this slows down the eval each time ?
          opts=utils.merge(opts,o);
      else
          opts=o;
      end
            
      N = b.N; M = numel(pts.x); k = b.k; eta=b.eta;
      ny = repmat(b.ny, [M 1]);
      d = repmat(pts.x, [1 N]) - repmat(b.y, [M 1]); % displacement matrix
      r = abs(d);

      if nargout==1,
          if eta==inf,
              % Value (no derivatives)
              A=utils.fundsol(r,k,'0',opts);
              if b.realflag, A=real(A); end;
              return              
          elseif eta==0,
              % Only first derivative
              [F0 A]=utils.fundsol(r,k,'1',opts);
              cosphi = -real(conj(ny).*d);
              A=A.*cosphi./r;
              if b.realflag, A=real(A); end;
              return
          else
              % Value and first derivative
              cosphi = -real(conj(ny).*d);
              [F0 F1]=utils.fundsol(r,k,'01',opts);
              A=F1.*cosphi./r+1i*eta*F0;
              if b.realflag, A=real(A); end;
              return
          end
      else
          rinv=1./r; dr=d.*rinv;
          if eta==inf,
              % Value and first derivative
              [A F1]=utils.fundsol(r,k,'01',opts);
              B=F1.*real(dr); C=F1.*imag(dr);
          elseif eta==0,
              % First and second derivative
              cosphi = -real(conj(ny).*d);
              [F0 F1 F2]=utils.fundsol(r,k,'12',opts);
              drr=dr.*rinv;
              A=F1.*rinv.*cosphi;
              B=F2.*real(drr).*cosphi-F1.*rinv.*(real(ny)+cosphi.*real(drr));
              C=F2.*imag(drr).*cosphi-F1.*rinv.*(imag(ny)+cosphi.*imag(drr));
          else
              % Value, first and second derivative
              cosphi = -real(conj(ny).*d); drr=dr.*rinv;
              [F0 F1 F2]=utils.fundsol(r,k,'012',opts);
              A=F1.*cosphi.*rinv+1i*eta*F0;
              B=F2.*real(drr).*cosphi-F1.*rinv.*(real(ny)+cosphi.*real(drr))+...
                  1i*eta*F1.*rinv.*real(d);
              C=F2.*imag(drr).*cosphi-F1.*rinv.*(imag(ny)+cosphi.*imag(drr))+...
                  1i*eta*F1.*rinv.*imag(d);
          end
          if b.realflag, A=real(A); B=real(B); C=real(C); end;
          if nargout==2,
              nx=repmat(pts.nx,[1 N]);
              nx=nx./abs(nx);
              B=real(nx).*B+imag(nx).*C;
          end
      end
      end
      
    function Nf = Nf(b)
      Nf = b.N;
    end
    
    function updateN(b,N)
    % UPDATEN - Update number N of MFS basis functions
        
        % Only possible if charge curve was given as function Z
        if isempty(b.Z), 
            warning('Cannot update number of MFS basis functions since charge curve Z is not given');
            return
        end        
        b.N=ceil(N*b.nmultiplier);
        % Create charge points, and if needed, outward-pointing unit normals...
        b.t=(1:b.N)/b.N;
        b.y=b.Z(b.t);
        if ~isempty(b.Zp), b.ny=-1i*b.Zp(b.t); b.ny = b.ny./abs(b.ny); end
    end
        

    function showgeom(bas, opts) % .................. crude show MFS pts, etc
    % SHOWGEOM - plots geometry of MFS basis charge points on current figure
    if nargin<2, opts = []; end
    h = ishold; hold on;                 % save hold state
      plot(real(bas.y), imag(bas.y), 'r+');
      if ~isempty(bas.ny)    % show dipole directions
        l = 0.1; ls = bas.y.'*[1 1]+l*bas.ny.'*[0 1]; % lines as separate rows
        plot(real(ls)', imag(ls)', 'm'); % flip to cols and plot them
      end
      if isfield(opts, 'label')
        n = ceil(numel(bas.y)/2);
        text(real(bas.y(n)), imag(bas.y(n)), opts.label);
      end
      if h, hold on; else, hold off; end  % restore hold state
    end % func
        
  end % methods
end
