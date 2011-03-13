% PROBLEM - abstract class defining interfaces for Helmholtz/Laplace BVPs/EVPs

% Copyright (C) 2008, 2009, Timo Betcke, Alex Barnett
classdef problem < handle
  properties
    segs                    % array of handles of segments in problem
    doms                    % array of handles of domains in problem
    k                       % overall wavenumber (or, if you will, frequency)
    A                       % BC inhomogeneity matrix (incl sqrt(w) quad wei)
    sqrtwei                 % row vec of sqrt of quadrature weights w
    bas                     % cell array of handles of basis sets in problem
    basnoff                 % dof index offsets for basis objs referred by bas
    N                       % total number of dofs in a problem
    co                      % basis coefficient (length N column vector)
  end
  
  methods % ------------------------------ methods common to problem classes
  
   function fillquadwei(pr)   % ............... fill quadrature weights vector
   % FILLQUADWEI - compute sqrt of bdry quadrature weights vector for a problem
   %
   % p.fillquadwei where p is a problem object fills the column vector p.sqrtwei
   %   with the square-root of quadrature weights for the segment BC / matching
   %   points in a problem. The ordering matches the row ordering of p.A matrix
   %   and p.rhs right-hand side vector.
   %
   % See also: BVP.SOLVECOEFFS
      pr.sqrtwei = [];        % get ready to stack stuff onto it as a big row
      for s=pr.segs
        if s.bcside==0          % matching condition
          w = sqrt(s.w);
          pr.sqrtwei = [pr.sqrtwei, w, w];  % stack twice the seg M
        elseif s.bcside==1 | s.bcside==-1   % BC (segment dofs natural order)
          pr.sqrtwei = [pr.sqrtwei, sqrt(s.w)]; 
        end
      end
      if isempty(pr.sqrtwei), warning('p.sqrtwei empty: probably you forgot to specify any BCs via setbc!'); end
    end % func
    
    function [N noff] = setupbasisdofs(pr)
    % SETUPBASISDOFS - set up indices of all basis degrees of freedom in problem
    %
    %   [N noff] = SETUPBASISDOFS(pr) sets problem basis set handle list pr.bas,
    %     overall # dofs pr.N, and problem basis object dof offsets pr.basnoff.
    %     It returns the last two. This is a helper routine for problem classes.
    %
    %    See also PROBLEM
      pr.bas = {};
      for d=pr.doms                  % loop over domains gathering basis sets
        pr.bas = {pr.bas{:} d.bas{:}};
      end
      pr.bas = utils.unique(pr.bas); % keep only unique basis handle objects
      noff = zeros(1, numel(pr.bas));
      N = 0;                         % N will be total # dofs (cols of A)
      for i=1:numel(pr.bas)
        noff(i) = N; N = N + pr.bas{i}.Nf; % set up dof offsets of bases
      end
      pr.N = N; pr.basnoff = noff;   % store stuff as problem properties
      if isempty(pr.bas), fprintf('warning: no basis sets in problem!\n'); end
    end
    
    function A = fillbcmatrix(pr, opts) %...........make bdry mismatch matrix
    % FILLBCMATRIX - computes matrix mapping basis coeffs to bdry mismatch
    %
    %   A = FILLBCMATRIX(pr) returns the matrix mapping all basis degrees of
    %   freedom in the problem to all segment boundary or matching condition
    %   inhomogeneity functions. With no output argument the answer is written
    %   to the problem's property A.
    %
    %   A = FILLBCMATRIX(pr, opts) allows certain options including:
    %   opts.doms: if present, restricts to contrib from only domain list doms.
    %   opts.trans: if present, gives translation applied to all segments before
    %    evaluation is done.
    %   opts.basobj: if present, basis objects (bas, basnoff) are read from
    %    problem or domain object basobj (which must have setupbasisdofs method)
    %    Note that in order to affect the problem segments they must be listed
    %    in some problem domain.
    %    (The above three are needed by blochmodeproblem and qpscatt)
    %
    %  Issues/Notes:
    %  * is now faster version, based on basis dofs rather than domain dofs
    %  * generalized to multiple affected domains per basis set, 6/4/10
      if nargin<2, opts=[]; end
      if ~isfield(opts, 'doms'), opts.doms = []; end
      if ~isfield(opts, 'trans'), trans = []; else trans = opts.trans; end
      if isempty(pr.sqrtwei), pr.fillquadwei; end
      if isfield(opts, 'basobj'), bob = opts.basobj; else bob = pr; end
      N = bob.setupbasisdofs;
      A = zeros(numel(pr.sqrtwei), N);       % A is zero when no basis influence
      m = 0;                                 % colloc index counter
      o = [];                                % opts to pass into basis.eval
      for s=pr.segs % ======== loop over segs
        % either use seg as target, or create a translated pointset (non-self):
        if isempty(trans), c = s; else c = pointset(s.x + trans, s.nx); end
        if s.bcside==0            % matching condition (2M segment dofs needed)
          ms = m + (1:2*size(s.x,1));     % 2M colloc indices for block row
          dp = s.dom{1}; dm = s.dom{2};   % domain handles on the + and - side
          for i=1:numel(bob.bas)           % all bases in problem or basis obj
            b = bob.bas{i};
            ns = bob.basnoff(i)+(1:b.Nf);     % dof indices for this bas
            talkp = utils.isin(b, dp.bas);   % true if this bas talks to + side
            talkm = utils.isin(b, dm.bas);   %                           - side
            if talkp ~= talkm                % talks to either one not both
              o.dom = dp; ind = 1;
              if talkm, o.dom = dm; ind = 2; end   % tell b.eval & a which side
              if isempty(opts.doms) | utils.isin(o.dom, opts.doms) % restrict?
                [Ab Abn] = b.eval(c, o);
                A(ms, ns) = repmat(pr.sqrtwei(ms).', [1 b.Nf]) .* [s.a(ind)*Ab; s.b(ind)*Abn]; % overwrite block, g stacked below f
              end
            elseif talkp & talkm             % bas talks to both (eg transm LP)
              o.dom = dp;
              if isempty(opts.doms) | utils.isin(o.dom, opts.doms) % restrict?
                [Ab Abn] = b.eval(c, o);  % + side contrib
                A(ms, ns) = A(ms, ns) + repmat(pr.sqrtwei(ms).', [1 b.Nf]) .* [s.a(1)*Ab; s.b(1)*Abn]; % add block, g stacked below f
              end
              o.dom = dm;
              if isempty(opts.doms) | utils.isin(o.dom, opts.doms) % restrict?
                [Ab Abn] = b.eval(c, o);  % - side contrib
                A(ms, ns) = A(ms, ns) + repmat(pr.sqrtwei(ms).', [1 b.Nf]) .* [s.a(2)*Ab; s.b(2)*Abn]; % add block
              end
            end
          end
            
        elseif s.bcside==1 | s.bcside==-1  % BC (M segment dofs, natural order)
          ind = (1-s.bcside)/2+1; % index 1,2 for which side the BC on (+,-)
          d = s.dom{ind};         % handle of domain on the revelant side
          ms = m+(1:size(s.x,1)); % colloc indices for this block row
          if isempty(opts.doms) | utils.isin(d, opts.doms) %restrict to domains?
            o = []; o.dom = d;      % b.eval may need to know in which domain
            for i=1:numel(bob.bas)   % all bases in problem or bas obj (not dom)
              b = bob.bas{i};      % (that's needed so i index reflects bob.bas)
              if utils.isin(b, d.bas)   % true if this bas talks to BC side
                if s.b==0               % only values needed, ie Dirichlet
                  Ablock = b.eval(c, o);
                  if s.a~=1.0, Ablock = s.a * Ablock; end
                else                    % Robin (includes Neumann)
                  [Ablock Anblock] = b.eval(c, o);
                  if s.a==0 & s.b==1.0  % Neumann
                    Ablock = Anblock;
                  else                  % Robin
                    Ablock = s.a*Ablock + s.b*Anblock;
                  end
                end
                A(ms,pr.basnoff(i)+(1:b.Nf)) = ...
                    repmat(pr.sqrtwei(ms).', [1 b.Nf]) .* Ablock; % write blk
              end
            end
          end
        end
        m = ms(end);                    % update counter in either case
      end % ================ loop over segs
      if nargout==0, pr.A = A; end     % this only stores internally if no outp
    end % func
    
    function [A Ax Ay] = evalbases(pr, p, opts) % .......... eval problem bases
    % EVALBASES - evaluate all basis sets in a domain object, on a pointset
    %
    %  A = EVALBASES(pr, p) returns matrix A whose jth column is the jth basis
    %   function in problem pr evaluated on the pointset p.
    %
    %  [A An] = EVALBASES(pr, p) returns also the normal derivatives using the
    %   normals associated with the pointset.
    %
    %  [A Ax Ay] = EVALBASES(pr, p) returns A and the basis x- and y-partial
    %   derivatives, ignoring the normals associated with the pointset
    %
    %  [A ...] = EVALBASES(pr, p, opts) allows options such as,
    %    opts.dom : overrides evaluation domain for bases (eg p = layerpot seg)
    %
    % Notes: 1) Although thought obsolete for domains, evalbases as a problem
    %   method is quite useful, esp. for periodic problems (nei>0, etc)
    %   2) For layer potential basis sets, and if p is also a segment
    %   object, jump relations will be taken into account if the opts.dom is
    %   specified.
    %   3) This routine can be viewed as PROBLEM.POINTSOLUTION without the
    %   final multiplication by coefficient vector.
    %   4) could make points lying outside all problem domains return NaN rows?
    %
    % See also: DOMAIN.EVALBASES, POINTSOLUTION, DOMAININDICES
      if nargin<3, opts = []; end
      if isempty(p.x), warning('evalbases on empty pointset!');
        A = []; Ax = []; Ay = []; return; end
      if nargout>1 && isempty(p.nx), error('pointset needs normals!'); end
      A = zeros(numel(p.x), pr.N);           % allocate
      if nargin>1, Ax = A; end; if nargin>2, Ay = A; end % only allocate if need
      di = NaN*zeros(size(p.x));             % NaN indicates in no domain
      if isfield(opts, 'dom'), doms = opts.dom; else doms = pr.doms; end
      for n=1:numel(doms)  % main loop is over domains, but bases are pr-indexed
        d = doms(n);
        ii = d.inside(p.x);
        if nnz(ii)==0, continue; end         % Do nothing if there's no elements
        di(ii) = n;                          % returns 1 if opts.dom override
        opts.dom = d;                        % b.eval might need know which dom
        pn = pointset(p.x(ii));              % points only in current dom
        if nargout>1, pn.nx = p.nx(ii); end
        for i=1:numel(pr.bas)                % loop over all bases...
          b = pr.bas{i}; ns = pr.basnoff(i)+(1:b.Nf);
          if utils.isin(b, d.bas)            % this bas talks to current dom?
            if nargout==1, A(ii,ns) = b.eval(pn, opts); % the 3 output styles
            elseif nargout==2, [Ad Adn] = b.eval(pn, opts);
              A(ii,ns) = Ad; Ax(ii,ns) = Adn;
            else, [Ad Adx Ady] = b.eval(pn, opts);
              A(ii,ns) = Ad; Ax(ii,ns) = Adx; Ay(ii,ns) = Ady; end
          end
        end
      end
    end
    
    function [u di] = pointsolution(pr, p, o) % ......eval soln on pointset
    % POINTSOLUTION - evaluate solution to a problem on a pointset, given coeffs
    %
    %  [u di] = pointsolution(pr, pts) returns array of values u, and
    %   optionally, domain index list di (integer array of same shape as u).
    %   Decisions about which domain a gridpoint is in are done using
    %   domain.inside, which may be approximate.
    %
    %  [u di] = pointsolution(pr, pts, opts) controls certain options, such as:
    %   opts.FMM = 1 uses Greengard-Gimbutas Helmholtz 2D FMM (wavenumber k>0),
    %                (certain bases only: mfsbasis, layerpots). DUMMY FOR NOW.
    %              0 uses direct filling of a dense matrix, hitting coeff vector
    %                (blocked according to opts.nmax, see Note 3, default 1e4)
    %
    %   Notes: 1) changed to reference dofs via bases rather than domains.
    %   2) A separate routine should be written for evaluation of u, u_n on
    %   boundary. For that matter, u_x and u_y should be accessible anytime too.
    %   3) splits up the computation if list of points is longer than opts.nmax.
    %   Ensures that not too much memory is eaten up by the computation.
    %   (Timo Betcke)
    %
    % See also GRIDSOLUTION.
      if nargin<3, o = []; end
      if ~isfield(o, 'nmax'), o.nmax = 3e3; end    % default blocking size
      % (since nmax=100 had big 20% speed hit due to oo-matlab overhead).
      if isfield(o, 'fmm'), o.FMM = o.fmm; end     % make case-insensitive
      if ~isfield(o, 'FMM'), o.FMM = 0; end        % default evaluation method
      if o.FMM, o.nmax = Inf; end                  % don't block the FMM!
      
      if isempty(pr.co), error('coefficient vector is empty!'); end
      if length(p.x)>o.nmax,
          Np=length(p.x);
          itcount=floor(Np/o.nmax);
          r=mod(Np, o.nmax);
          u=zeros(Np,1); di=zeros(Np,1);
          for j=1:itcount+1,
              if j<=itcount,
                  indrange=(j-1)*o.nmax+1:j*o.nmax;
              elseif r>0,
                  indrange=(j-1)*o.nmax:Np;
              else break
              end
              if ~isempty(p.nx),             
                  ptemp=pointset(p.x(indrange),p.nx(indrange));
              else
                  ptemp=pointset(p.x(indrange));
              end
              [ut,dit]=pr.pointsolution@problem(ptemp, o); % recursive!(depth 1)
              u(indrange)=ut; di(indrange)=dit;
          end
          return
      end
      di = NaN*zeros(size(p.x));                    % NaN indicates in no domain
      u = di;                                       % solution field
      for n=1:numel(pr.doms)
        d = pr.doms(n);
        ii = d.inside(p.x);
        if nnz(ii)==0, 
            continue
        end                                % Do nothing if there are no elements
        di(ii) = n;
        u(ii) = 0;                           % accumulate contribs to u in dom
        opts.dom = d;                        % b.eval might need know which dom
        for i=1:numel(pr.bas)                % loop over all bases...
          b = pr.bas{i};
          if utils.isin(b, d.bas)            % this bas talks to current dom?
            co = pr.co(pr.basnoff(i)+(1:b.Nf)); % extract coeff vec for basis
            if o.FMM & b.HFMMable & b.k>0       % Helmholtz only, not k=0!
              u(ii) = u(ii) + b.evalFMM(co, p.x(ii));
            else
              Ad = b.eval(pointset(p.x(ii)), opts); % fill dense matrix
              u(ii) = u(ii) + Ad * co;           % add effect of this basis obj
            end
          end
        end
      end
    end % func

    function di = domainindices(pr, p)
    % DOMAININDICES - return domain (problem) indices of points in pointset
    %
    % Note: just a helper routine extracting the di output part of pointsolution
      di = NaN*zeros(size(p.x));              % NaN indicates in no domain
      for n=1:numel(pr.doms)
        d = pr.doms(n);
        ii = d.inside(p.x);
        di(ii) = n;
      end
    end
    
    function [u gx gy di] = gridsolution(pr, o) % ......... eval soln on grid
    % GRIDSOLUTION - evaluate solution to a problem over a grid, given coeffs
    %
    %  [u gx gy di] = gridsolution(pr, opts) returns array of values u, and
    %   optionally, x- and y-grids (1D lists) gx, gy, and domain index list di
    %   (integer array of same shape as u). Decisions about which domain a
    %   gridpoint is in are done using domain.inside, which may be approximate.
    %
    %  Other options in opts are passed to pointsolution; see its doc page.
    %
    % To do: * keep evalbases matrices Ad for later use, multiple RHS's etc.
    % * what if evalbases matrices too big to store, sum basis vals by hand?
    %
    % See also POINTSOLUTION, GRIDBOUNDINGBOX
      o = pr.gridboundingbox(o);
      gx = o.bb(1):o.dx:o.bb(2); gy = o.bb(3):o.dx:o.bb(4);  % plotting region
      [xx yy] = meshgrid(gx, gy); zz = xx + 1i*yy;  % keep zz rect array
      [u di] = pr.pointsolution(pointset(zz(:)), o); % make zz a col vec
      u = reshape(u, size(xx));
      di = reshape(di, size(xx));
    end % func
    
    function o = gridboundingbox(pr, o)
    % GRIDBOUNDINGBOX - set default options giving grid spacing and bound box
    %
    %  This code was pulled from gridsolution so scattering class can access it
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
      tiny = 1e-12;                      % infinitesimal jog for stability
      o.bb(1) = o.bb(1) + tiny; o.bb(3) = o.bb(3) + tiny; % jog grid (inside) 
    end
    
    function h = showbdry(pr, o)   % ........................ crude plot bdry
    % SHOWBDRY - shows boundary segments in a problem, with their natural sense
      if nargin<2, o = []; end
      h = domain.showsegments(pr.segs, ones(size(pr.segs)), o);
    end
    
    function h = showbasesgeom(pr)   % ................ crude plot bases geom
    % SHOWBASESGEOM - plot geometry of all basis objects in a problem object
      h = [];   % dummy graphics handle for now
      pr.setupbasisdofs;
      for i=1:numel(pr.bas)
        opts.label = sprintf('%d', i);              % label by problem's bas #
        pr.bas{i}.showgeom(opts);
      end
    end
    
    function h = plot(pr)
    % PLOT - show all geometry in a problem object: segments, basis geoms...
      h = [pr.showbdry; pr.showbasesgeom];
    end
    
    function setoverallwavenumber(pr, k) % ................. overall k
    % SETOVERALLWAVENUMBER - set problem k in each domain using refractive index
    %
    %  setoverallwavenumber(pr, k) propagates overall wavenumber k to each
    %   domain, and its basis sets. If a domain has refractive index n, then
    %   its wavenumber will become nk.
      if isnan(k), error('k must be a number!'); end
      pr.k = k;
      for d=pr.doms
        if isnan(d.refr_ind), error('each domain index must be a number!'); end
        d.k = d.refr_ind * k;
      end
    end
    
    function updateN(pr,N)
    % UPDATEN - Update the number of basis functions in each basis set
    %
    % updateN(pr,N) sets the degree/number of basis functions in each
    %   basis set of the problem pr to be N times the basis set's nmultiplier
    %   property. This is useful for convergence studies when all basis degrees
    %   should be changed in proportion.
        for d=pr.doms,
            for j=1:length(d.bas),
                d.bas{j}.updateN(N); 
            end
        end
        pr.setupbasisdofs;       % refreshes pr.N and pr.basnoff
    end
    
    
    function [u gx gy di h] = showsolution(pr, o) % ............. plot solution
    % SHOWSOLUTION - plot figure with solution Re u over all domains in problem
    %
    % showsolution(pr) plots an image of the solution field u (its real part
    %   if complex) over a rectangular grid covering the domains in the
    %   problem pr.
    %
    % showsolution(pr, opts) also passes in an option structure containing
    %   optional fields:
    %      opts.imag = true, plots imag instead of real part
    %      opts.bdry = true, shows boundary too (default)
    %      opts.comparefunc = handle to function to subtract pointwise
    %                         (this function must cope with vector inputs)
    %   Other fields are passed to problem.gridsolution (see its documentation)
    %
    % [u gx gy di h] = showsolution(pr, ...) also returns solution data in
    %   same format at problem.gridsolution, and h is a handle to the image
    %   graphical object
    %
    %  Issues: * make a real/complex flag
    %          * store Ad eval matrices for later access (expensive to fill)?
    %          * Make opts.comparefunc, if a vector of func handles, apply
    %            different handle to each domain in the problem.
    %
    % Also see: GRIDSOLUTION, SCATTERING.SHOWTHREEFIELDS
      if nargin<2, o = []; end
      if ~isfield(o, 'imag'), o.imag = 0; end
      if ~isfield(o, 'bdry'), o.bdry = 0; end
      
      [u gx gy di] = pr.gridsolution(o);   % expensive: do basis evaluations
      
      if isfield(o, 'comparefunc')
        if isa(o.comparefunc, 'function_handle')
          [xx yy] = meshgrid(gx, gy); zz = xx + 1i*yy;
          u = u - o.comparefunc(zz);       % eval comparison func over grid
        else
          error('opts.comparefunc must be a function_handle!');
        end
      end
      % note use of transparency outside all domains...
      if o.imag, h = imagesc(gx, gy, imag(u), 'alphadata', ~isnan(di));
        %title('Im[u]');
      else, h = imagesc(gx, gy, real(u), 'alphadata', ~isnan(di));
        %title('Re[u]');
      end
      utils.goodcaxis(u); axis equal tight; colorbar; colormap(jet(256));
      set(gca,'ydir','normal'); hold on;
      if o.bdry, pr.showbdry; end
    end % func
    
    % *** Methods to be written ........... ****
    [u un] = bdrysolution(pr, seg, pm) % ........... evaluate soln on a bdry
    
  end % methods
   
  methods (Abstract) % ------------------------------------------------- 
    % none
  end
end
