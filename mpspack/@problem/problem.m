% PROBLEM - abstract class defining interfaces for
%  Helmholtz or Laplace BVP or EVP.

classdef problem < handle
  properties
    segs                    % handles of segments in BVP
    doms                    % handles of domains in BVP
    k                       % overall wavenumber
    A                       % BC inhomogeneity matrix (incl sqrt(w) quad wei)
    sqrtwei                 % row vec of sqrt of quadrature weights w
    co                      % basis coefficients
  end
  
  methods % ------------------------------ methods common to problem classes
  
   function fillquadwei(pr)   % ............... fill quadrature weights vector
   % FILLQUADWEI - compute sqrt of bdry quadrature weights vector for a problem
      pr.sqrtwei = [];        % get ready to stack stuff onto it as a big row
      for s=pr.segs
        if s.bcside==0          % matching condition
          w = sqrt(s.w);
          pr.sqrtwei = [pr.sqrtwei, w, w];  % stack twice the seg M
        elseif s.bcside==1 | s.bcside==-1   % BC (segment dofs natural order)
          pr.sqrtwei = [pr.sqrtwei, sqrt(s.w)]; 
        end
      end
    end % func
    
    function A = fillbcmatrix(pr)   % ........... make bdry discrepancy matrix
    % FILLBCMATRIX - computes matrix mapping basis coeffs to bdry discrepancy 
      if isempty(pr.sqrtwei), error('must fill quadrature weights first'); end
      N = 0;                        % N will be total # dofs (cols of A)
      for d=pr.doms, d.noff = N; N = N+d.Nf; end % setup noff's in d's
      A = [];            % get ready to stack block rows, dof order matches rhs
      for s=pr.segs
        if s.bcside==0            % matching condition (2M segment dofs needed)
          % to do...
          
          
        elseif s.bcside==1 | s.bcside==-1  % BC (M segment dofs, natural order)
          ind = (1-s.bcside)/2+1; % index 1 or 2 for which side the BC on
          d = s.dom(ind);         % handle of domain on the revelant side
          d = d{1};               % ugly, extracts domain handle from cell
          if s.b==0               % only values needed, ie Dirichlet
            Ablock = s.a * d.evalbases(s);
          else;                   % Robin (includes Neumann)
            [Ablock Anblock] = d.evalbases(s);
            Ablock = s.a*Ablock + s.b*Anblock;
          end
          Arow = zeros(size(s.x, 1), N);     % start with empty block-row
          %size(Arow), size(Ablock), d.noff, d.Nf   % debug
          Arow(:,d.noff+(1:d.Nf)) = Ablock;  % copy in nonzero block
          A = [A; Arow];                     % stack block rows
        end
      end
      A = A .* repmat(pr.sqrtwei.', [1 N]); % quad wei: better do this in bits?
      if nargout==0, pr.A = A; end     % this only stores internally if no outp
    end % func
    
    function [u di] = pointsolution(pr, p) % ...........eval soln on pointset
    % POINTSOLUTION - evaluate solution to a problem on a pointset, given coeffs
    %
    %  [u di] = gridsolution(pr, pts) returns array of values u, and
    %   optionally, domain index list di (integer array of same shape as u).
    %   A separate routine should be used for evaluation of u, u_n on boundary.
    %   Decisions about which domain a gridpoint is in are done using
    %   domain.inside, which may be approximate.
    %
    % See also GRIDSOLUTION.
      di = NaN*zeros(size(p.x));                    % NaN indicates in no domain
      u = di;                                       % solution field
      for n=1:numel(pr.doms)
        d = pr.doms(n);
        ii = d.inside(p.x);
        di(ii) = n;
        Ad = d.evalbases(pointset(p.x(ii))); % bases eval matrix for pts in dom
        co = pr.co(d.noff+(1:d.Nf));         % extract coeff vector for domain d
        u(ii) = Ad * co;                     % set only values inside
      end
    end % func
    
    function [u gx gy di] = gridsolution(pr, o) % ......... eval soln on grid
    % GRIDSOLUTION - evaluate solution to a problem over a grid, given coeffs
    %
    %  [u gx gy di] = gridsolution(pr, opts) returns array of values u, and
    %   optionally, x- and y-grids (1d lists) gx, gy, and domain index list di
    %   (integer array of same shape as u). Decisions about which domain a
    %   gridpoint is in are done using domain.inside, which may be approximate.
    %
    % To do: * keep evalbases matrices Ad for later use, multiple RHS's etc.
    % * what if evalbases matrices too big to store, sum basis vals by hand?
    %
    % See also POINTSOLUTION.
      if nargin<2, o = []; end
      if ~isfield(o, 'dx'), o.dx = 0.03; end    % default grid spacing
      if o.dx<=0, error('dx must be positive!'); end
      if ~isfield(o, 'bb')                      % default bounding box
        bb = [];
        for d=pr.doms
          bb = [bb; d.boundingbox];
        end
        o.bb([1 3]) = min(bb(:,[1 3]), [], 1);  % find box enclosing all BBs
        o.bb([2 4]) = max(bb(:,[2 4]), [], 1);
      end
      o.bb(1) = o.dx * floor(o.bb(1)/o.dx);         % quantize to grid through 0
      o.bb(3) = o.dx * floor(o.bb(3)/o.dx);         % ... make this optional?
      gx = o.bb(1):o.dx:o.bb(2); gy = o.bb(3):o.dx:o.bb(4);  % plotting region
      [xx yy] = meshgrid(gx, gy); zz = xx + 1i*yy;  % keep zz rect array
      [u di] = pr.pointsolution(pointset(zz));
    end % func
    
    function h = showbdry(pr)   % ........................ crude plot bdry
    % SHOWBDRY - shows boundary segments in a problem, with their natural sense
      h = domain.showsegments(pr.segs, ones(size(pr.segs)));
    end
    
    
    % *** Methods to be written ........... ****
    [u un] = bdrysolution(pr, seg, pm) % ........... evaluate soln on a bdry
    
    
  end % methods
   
  methods (Abstract) % ------------------------------------------------- 
    % none
  end
end
