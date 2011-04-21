% BVP - define a Helmholtz or Laplace boundary value problem from domains
%
%  pr = bvp(doms) creates a BVP object pr which comprises the domains doms
%   (doms is a row vec of domain handles), using their segments and the
%   boundary conditions they carry, and the basis sets associated with the
%   domains.

% Copyright (C) Alex Barnett & Timo Betcke 2008-2011.

classdef bvp < problem & handle
  properties
    rhs                             % right hand side (sqrt(w).f) col vector
  end
  
  methods % ---------------------------------- methods particular to BVPs
    function pr = bvp(doms) % .................... constructor
      if nargin==0, return; end    % needs empty constructor
      pr.doms = doms;
      pr.segs = [];             % now build a list of segments, crude O(N^2)
      for d=doms
        for s=d.seg
          if isempty(find(pr.segs==s))  % only add the segment if it's new
            pr.segs = [pr.segs s];
          end
        end
      end
    end

    function rhs = fillrighthandside(pr)   % ............... make RHS vector
      if isempty(pr.sqrtwei), error('must fill quadrature weights first'); end
      rhs = [];                 % get ready to stack stuff onto it as a big col
      for s=pr.segs
        if s.bcside==0          % matching condition
          if isnumeric(s.f)
            rhs = [rhs; s.f; s.g];           % data vec, stack as one big col
          else
            rhs = [rhs; s.f(s.t); s.g(s.t)]; % func of t, stack as one big col
          end
        elseif s.bcside==1 | s.bcside==-1    % BC (segment dofs natural order)
          if isnumeric(s.f)
            rhs = [rhs; s.f];       % data vector, stack as one big column
          else
            rhs = [rhs; s.f(s.t)];  % function of t, stack as one big column
          end
        end
      end
      rhs  = rhs .* pr.sqrtwei.';          % include sqrt quad weights factor
      if nargout==0, pr.rhs = rhs; end     % this only stores pr.rhs if no outp
    end
    
    function co = linsolve(pr, o) % ........................ linear solve
    % LINSOLVE - solve the full linear system matrix directly or iteratively
    %
    % pr.linsolve where pr is a problem (bvp, etc) solves the linear system
    %  Ax = b where b is the right hand side (pr.rhs) and A the system matrix
    %  (pr.A, or this may instead be a function that maps x to Ax). Solution
    %  x is stored as the property pr.co (solution coefficients).
    %
    % co = pr.linsolve outputs instead of storing the the problem property.
    %
    % pr.linsolve(opts) controls optional settings:
    %    opts.meth = 'direct' or 'iter' controls solution method: dense
    %                factorization, or GMRES iterations.
    %    opts.eps (default 1e-12) controls desired residual error in GMRES
    %    opts.matname (default 'A') is text string giving the name of the
    %                system matrix (a problem property).
    %
    % Notes:
    % 1) move this up to @problem ?
      if nargin<2, o = []; end
      if ~isfield(o, 'meth'), o.meth = 'direct'; end    % default method
      if ~isfield(o, 'eps'), o.eps = 1e-12; end         % desired residual
      if ~isfield(o, 'matname'), o.matname = 'A'; end   % point to pr.A
      A = eval(['pr.' o.matname]);                      % system matrix/applier
      
      if isempty(A), error('system matrix is empty!'); end
      
      if strcmp(o.meth, 'direct')
        if ~isnumeric(A)
          error('cannot direct solve given only an operator applier function!');
        end
        disp('direct solve...')
        co = A \ pr.rhs;
      elseif strcmp(o.meth, 'iter') % NB here o.mat can be matrix or func...
        maxit = min(1e3, numel(pr.rhs));
        disp('GMRES solve...')
        [co] = gmres(A, pr.rhs, [], o.eps, maxit);
        %[co flag relres iter] = gmres(pr.A,pr.rhs,[],o.eps,maxit); % no verb
        %fprintf('iteration numbers: %d (inner %d)\n', iter(1), iter(2))
        % do some error handling!
        % ...
      else, error('unknown meth');
      end
      if nargout==0, pr.co = co; end   % only store internally if no output
    end
    
    function co = solvecoeffs(pr, o) % ................. fill and linear solve
    % SOLVECOEFFS - solve a BVP for basis coefficients (possibly after filling)
    %
    %  co = solvecoeffs(pr) where pr is a BVP object, sets up the RHS and the
    %   system matrix, then does a dense linear solve (which is the default).
    %
    %  co = solvecoeffs(pr, opts) lets the user control the solution options:
    %    opts.FMM = 0 or 1 controls if a dense matrix or FMM is used
    %                 to apply the operator
    %    opts.meth = 'factor' or 'iter' controls solution method, dense
    %                factorization, or GMRES iterations.
    %
    %  See also: FILLBCMATRIX, APPLYBCMATRIX, LINSOLVE
    
      if nargin<2, o = []; end
      if ~isfield(o, 'FMM'), o.FMM = 0; end             % default applier
      
      if isempty(pr.sqrtwei), pr.fillquadwei; end
      if isempty(pr.rhs), pr.fillrighthandside; end
      pr.setupbasisdofs;
      
      % set up operator representation type...
      if o.FMM==0
        pr.fillbcmatrix;   % fills pr.A with dense matrix
      else
        pr.A = @(co) pr.applybcmatrixFMM(co, o); % operator applier function
        if ~isfield(o, 'meth'), o.meth = 'iter'; end    % default for FMM
      end
      co = pr.linsolve(o); % do the linear solve (o controls which type)
      if nargout==0, pr.co = co; end   % only store internally if no output
    end
    
    function r = bcresidualnorm(pr, o)    % ............... L2 bdry error
    % BCRESIDUALNORM - norm of mismatch residual for coefficient vector in a BVP
      if nargin<2, o = []; end
      if isnumeric(pr.A)
        y = pr.A * pr.co;
      else
        y = pr.applybcmatrixFMM(pr.co, o);
      end
      r = norm(y - pr.rhs);
    end
    
  end
end
