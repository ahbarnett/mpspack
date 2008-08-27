% BLOCHMODEPROBLEM - problem class of finding quasi-periodic (Bloch) eigenmodes
%
%   pr = BLOCHMODEPROBLEM(uc, cond, d) makes a problem object with unit cell uc,
%    a single connected domain cond, and inclusion (disconnected) domains d
%    (may be an array of domain handles).
%
%   See also QPUNITCELL

classdef blochmodeproblem < bvp & problem
  properties
    uc                  % unit cell for principal connected `air' component
    sourcenei           % controls # src neighbors: 0 gives 1x1, 1 gives 3x3 ...
    mismatchnei         % controls # mismatch mixture neighbors (same as above)
    Adata               % struct giving BC matrix for local copies, etc
    Q                   % quasi-periodizing basis to unit cell discrep matrix
    C                   % obstacle basis to unit cell discrep matrix
    B                   % quasi-periodizing basis to obstacle mismatch matrix
    wco                 % whole-problem coeff vector
  end
  methods
    function pr = blochmodeproblem(uc, condom, doms, opts) % ...... constructor
      pr = pr@bvp([condom doms]);      % make an instance of a BVP (see Docu)
      % note that pr.doms(1) is always the connected domain (condom)
      if ~isempty([condom doms])
        if numel(condom)~=1
          error('must be exactly one lattice-connected domain!')
        end
        condom.isair = 1;            % flag the connected domain as `air'
        for d=doms, d.isair = 0; end
      end
      pr.uc = uc;
      pr.k = uc.k;
      if isempty(uc.bas)
        fprintf('warning: qpunitcell uc has no bases!\n');
      end
      if nargin<4, opts = []; end
      if ~isfield(opts, 'sourcenei') | isempty(opts.sourcenei)
        opts.sourcenei = 0; end                                 % default
      if ~isfield(opts, 'mismatchnei') | isempty(opts.mismatchnei)
        opts.mismatchnei = 0; end                               % default
      pr.sourcenei = opts.sourcenei; pr.mismatchnei = opts.mismatchnei;
    end
    
    function h = showbdry(pr) % ...................... overloads problem
    % SHOWBDRY - shows boundary segments and unit cell in Bloch mode problem
      h = domain.showsegments(pr.segs, ones(size(pr.segs)));
      h = [h; pr.uc.plot];
    end
    
    function C = fillobstodiscrep(pr) % ......................... C
    % obstacle basis to unit-cell discrepancy matrix
      pr.setupbasisdofs;            % note, only for obstacle + condom
      uc = pr.uc;
      M = numel(vertcat(uc.seg.x));
      C = zeros(M, pr.N);
      nei = pr.sourcenei;
      opts.dom = pr.doms(1);               % eval in the condom for obstacle
      for i=1:numel(pr.bas)
        b = pr.bas{i};
        if utils.isin(b, pr.doms(1).bas)  % does this b bas affect conn. dom? 
          ns = pr.basnoff(i)+(1:b.Nf);    % dof indices for this bas
          opts.nei = nei;                 % tell eval how many neighbor copies
          C(:,ns) = b.evalunitcelldiscrep(uc, opts);
          %C = C .* repmat(pr.uc.discrepsqrtwei.', [1 pr.N]);% no L2 discrep wei
        end
      end
      if nargout==0, pr.C = C; end     % this only stores internally if no outp
    end % func
    
    function A = fillbcmatrix(pr) % ......................... overloads problem
    % FILLBCMATRIX - fills A matrix including phased local src copies, stores
    %
    %   For the phased neighbor copies, only considers effect of conn. dom
    %   bases on segments touching conn. dom. Stores unphased copy matrices
    %   and reuses (with new phases) if considers they are still current.
      nei = pr.sourcenei; uc = pr.uc;
      wantcopies = isempty(pr.Adata);   % need to compute A or its copies?
      if ~wantcopies                    % various reasons to need to recompute:
        wantcopies = (pr.k~=pr.Adata.k | pr.N~=size(pr.Adata.A, 2) | ...
                      (2*pr.sourcenei+1)^2~=size(pr.Adata.A, 3));
      end
      if isempty(pr.sqrtwei), pr.fillquadwei; end
      N = pr.setupbasisdofs;
      if wantcopies
        pr.Adata.k = pr.k;
        pr.Adata.A = zeros([numel(pr.sqrtwei) N (2*nei+1)^2]);  % alloc copies
      end
      fprintf('\t\twantcopies = %d, nei = %d\n', wantcopies, nei)
      o = []; o.doms = pr.doms(1);  % opts to only eval contrib from conn. dom
      c = 1;                   % counter of copies
      A = zeros([numel(pr.sqrtwei) N]);    % always zero the rephased sum matrix
      for n=-nei:nei
        for m=-nei:nei
          if wantcopies
            if n==0 & m==0     % the original, no restriction on domains
              pr.Adata.A(:,:,c) = pr.fillbcmatrix@bvp; % superclass meth no opts
            else               % one of local neighbor copies, use only condom
              o.trans = -n*uc.e1 - m*uc.e2; % translate target from original
              pr.Adata.A(:,:,c) = pr.fillbcmatrix@bvp(o);
            end
          end
          A = A + uc.a^n * uc.b^m * pr.Adata.A(:,:,c);
          c = c+1;
        end
      end
      if nargout==0, pr.A = A; end     % this only stores internally if no outp
    end
    
  %  function B = fillqptomismatch(pr) % ......................... B
    % FILLQPTOMISMATCH - quasi-periodizing basis to obstacle mismatch matrix B
   %   N = uc.setupbasisdofs;
   %   B = zeros([numel(pr.sqrtwei) N]);    % always zero the rephased sum matrix
   %   for n=-nei:nei
   %     for m=-nei:nei
   %       o.trans = 
   %       pr.fillbcmatrix(o);
   % 
   %     end
   %   end
   %   
   % end
    
    function B = fillqptomismatch(pr) % ......................... B
    % quasi-periodizing basis to obstacle mismatch matrix - lots of code
    % in common with problem.fillbcmatrix, could break out into functions?
      if isempty(pr.sqrtwei), pr.fillquadwei; end
      uc = pr.uc;
      [N noff] = uc.setupbasisdofs;
      B = zeros(numel(pr.sqrtwei), N);
      nei = pr.mismatchnei;
      m = 0;
      for s=pr.segs
        if s.bcside==0                   % matching
          ms = m + (1:2*size(s.x,1));     % 2M colloc indices for block row
          % ...
          
        elseif s.bcside==1 | s.bcside==-1  % BC (M segment dofs, natural order)
          ind = (1-s.bcside)/2+1; % index 1,2 for which side the BC on (+,-)
          d = s.dom{ind};         % handle of domain on the revelant side
          o = []; o.dom = d;      % b.eval may need to know in which domain
          ms = m+(1:size(s.x,1)); % M colloc indices for this block row
          if d==pr.doms(1)        % only do if is a BC of the connected domain
            for i=1:numel(uc.bas)        % all QP bases in unit cell
              b = uc.bas{i};
              ns = noff(i)+(1:b.Nf);     % dof indices for this bas
              if nei==0
                if s.b==0               % only values needed, ie Dirichlet
                  Bb = b.eval(s, o);
                  if s.a~=1.0, Bb = s.a * Bb; end
                else                    % Robin (includes Neumann)
                  [Bb Bnb] = b.eval(s, o);
                  if s.a==0 & s.b==1.0  % Neumann
                    Bb = Bnb;
                  else                  % Robin
                    Bb = s.a*Bb + s.b*Bnb;
                  end
                end
                B(ms, ns) = repmat(pr.sqrtwei(ms).', [1 b.Nf]) .* Bb;% write blk
              else                      % neighborhood mixture of mismatch
                % ...
              end
            end
          end
        end
        m = ms(end);                    % update counter in either case        
      end % segs loop
      if nargout==0, pr.B = B; end     % this only stores internally if no outp
    end
    
    function M = wholematrix(pr, opts)
    % WHOLEMATRIX - compute full obstacle-and-quasiperiodizing matrix
      if nargin<2, opts=[]; end
      if ~isfield(opts, 'twonorm'), opts.twonorm = 0; end % true if L2 weights
      pr.fillquadwei; pr.sqrtwei = ones(size(pr.sqrtwei)); % vital, no left mult
      pr.fillbcmatrix;
      % fac of 2 on first block row takes I/2 to I...
      M = [2*pr.A, 2*pr.fillqptomismatch; ...
           pr.fillobstodiscrep, pr.uc.evalbasesdiscrep];
      if opts.twonorm
        M = pr.layerpottwonormreweight(M);
      end     
      if nargout==0, pr.M = M; end     % this only stores internally if no outp
   end

    function [u gx gy di] = gridsolution(pr, o) % ......... eval soln on uc grid
    if ~isfield(o, 'dx'), o.dx = 0.03; end    % default grid spacing
      if o.dx<=0, error('dx must be positive!'); end
      if ~isfield(o, 'bb')                      % default bounding box
        o.bb = pr.uc.boundingbox;
      end
      o.bb(1) = o.dx * floor(o.bb(1)/o.dx);         % quantize to grid through 0
      o.bb(3) = o.dx * floor(o.bb(3)/o.dx);         % ... make this optional?
      gx = o.bb(1):o.dx:o.bb(2); gy = o.bb(3):o.dx:o.bb(4);  % plotting region
      [xx yy] = meshgrid(gx, gy); zz = xx + 1i*yy;  % keep zz rect array
      ii = pr.uc.inside(zz);
      u = NaN*zeros(size(zz)); di = u;         % stuff outside the uc
      [u(ii) di(ii)] = pr.pointsolution(pointset(zz(ii)));
    end
    
    function [u di] = pointsolution(pr, p) % ........... overloads from problem
    %  note that includes rephased multiple source copies (done by moving targ)
    %
      pr.co = pr.wco(1:pr.N);                % obstacle dofs coeff vec
      nei = pr.sourcenei; uc = pr.uc;
      u = zeros(size(p.x));
      for n=-nei:nei
        for m=-nei:nei
          pc = pointset(p.x -n*uc.e1 - m*uc.e2, p.nx);  % fake by moving target
          [unm di] = pointsolution@problem(pr, pc);     % call superclass method
          u = u + uc.a^n * uc.b^m * unm;                % source phases
        end
      end
      qpco = pr.wco(pr.N+1:end);             % QP dofs coeff vec
      d = pr.doms(1);                        % the condom
      ii = d.inside(p.x);
      opts.dom = d;                          % b.eval might need know which dom
      for i=1:numel(pr.uc.bas)               % loop over all QP bases...
        b = pr.uc.bas{i};
        Ad = b.eval(pointset(p.x(ii)), opts);
        co = qpco(pr.uc.basnoff(i)+(1:b.Nf)); % extract coeff vector for basis
        u(ii) = u(ii) + Ad * co;              % add contrib from this basis obj
      end
    end % func

    function M = layerpottwonormreweight(pr, M) % ..... reweights LP op M for L2
      if nargin<2, M = pr.M; end
      pr.fillquadwei; 
      sqrtw = [pr.sqrtwei, pr.uc.discrepsqrtwei]; % row vector
      M = M .* repmat(1./sqrtw, [size(M,1) 1]); % kills a sqrtw factor on right
      M = M .* repmat(sqrtw.', [1 size(M,2)]); % new sqrtw factor on left
      if nargout==0, pr.M = M; end     % this only stores internally if no outp
    end % func
    
    function wholecoeffvectwonormreweight(pr)
      pr.fillquadwei; 
      sqrtw = [pr.sqrtwei, pr.uc.discrepsqrtwei]; % row vector
      pr.wco  = pr.wco ./ sqrtw.';           % returns to density func values
    end % func
  end % methods
end
