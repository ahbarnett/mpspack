% EVP - define a Laplace eigenvalue problem for a single interior domain
%
% p = evp(dom) creates an EVP object from the domain dom, comprising the
%   domain along with its segments the boundary conditions they carry, and
%   the basis sets associated with the domain.
%
%   Preliminary version: only one domain is supported, with Dirichlet BCs,
%   and the only solution method is 
%   root-search on the Fredholm determinant of Id minus double-layer operator.

% Copyright (C) 2010 Alex Barnett, Timo Betcke

classdef evp < problem & handle
  properties
    kj                        % eigenwavenumbers = sqrt(eigenvalues) (col vec)
    coj                       % basis coeffs (columns of array) of modes
    err                       % struct for error estimates errors of kj:
                              %   err.ej - rootfinding errors (col vec)
                              %   err.minsigj - min op sing vals (col vec)
    ndj                       % normal derivatives (cols of array) of modes
    kwin                      % wavenumber window in which kj were requested
  end

  methods % ----------------------------------------- particular to EVP
    function pr = evp(doms) % .........................constructor
      if nargin==0, return; end    % needs empty constructor
      if numel(doms)>1, error('only single domain currently supported!'); end
      pr.doms = doms; % Following code duplicates bvp, should be farmed out:
      pr.segs = [];             % now build a list of segments, crude O(N^2)
      for d=doms
        for s=d.seg
          if isempty(find(pr.segs==s))  % only add the segment if it's new
            pr.segs = [pr.segs s];
          end
        end
      end
    end
  
    function [kj err coj ndj] = solvespectrum(p, kwin, meth, o)
    % SOLVESPECTRUM - compute spectrum, and possibly modes, in wavenumber window
    %
    % [kj] = SOLVESPECTRUM(p, [klo khi]) uses a default method to find all
    %   eigenwavenumbers kj lying in the wavenumber interval [klo, khi], for the
    %   eigenvalue problem object p.
    %   Laplacian eigenvalues are simply the squares of the kj.
    %   The problem also stores the kj list as p.kj
    %
    % [kj err coj ndj] = ... also returns basis coeffcients, normal derivatives,
    %   and error estimates struct, if possible.
    %
    % [...] = SOLVESPECTRUM(p, [klo khi], meth, opts) allows one to choose the
    %   method,
    %     meth = 'fd' : Fredholm det (layer-potential bases only)
    %                   Error estimates are dist to real axis of Boyd roots
    %   and choose options,
    %     opts.verb = 0,1,... : verbosity (0=silent, 1=diagnostic text, ...)
    %     opts.modes = 0,1 : if true, compute modes (storing in p)
    %
    % Currently the only solution method is root-search on the Fredholm
    %  determinant of Id minus double-layer operator, for a single domain, with
    %  Dirichlet BCs everywhere on the boundary.
    %
    % Notes / issues:
    % * Still very alpha code!
    % * keep more data of all the det evals, pass out?
      if nargin<3 | isempty(meth), meth='fd'; end   % default method
      if nargin<4, o = []; end
      if ~isfield(o, 'tol'), o.tol = 1e-8; end   % requested tolerance
      if ~isfield(o, 'verb'), o.verb = 1; end; v = o.verb; % default verbosity
      if ~isfield(o, 'modes'), o.modes = 0; end
      if numel(p.doms)~=1, error('problem must have exactly one domain!'); end
      d = p.doms;                            % the domain
      wantmodes = nargout>2 | o.modes;
      klo = kwin(1); khi = kwin(2); p.kwin = kwin;
      p.fillquadwei; M = numel(p.sqrtwei);   % total # boundary pts
      if v, fprintf('mean M ppw @ khi = %.3g\n', 2*pi*M/khi/d.perim); end
      
      if strcmp(meth,'fd') % ............... Fredholm det method
        k = klo; kj = []; ej = []; N = M;  % square system
        kscale = 70/d.diam;                % est k at which level density ~ 20
        io.Ftol = 1e-10; io.tol = o.tol;   % allows leeway in collecting roots
        io.disp = v;                       % opts for interval search
        while k<khi
          ktop = min(k + kscale*min(1/k,0.2), khi); % top of current k window
          % option to reset N here, requadrature all segs (in proportion?)
          % ...
          %f = @(k) det(eye(N) - 2*layerpot.D(k, s)); % for one closed segment
          f = @(k) det(p.fillfredholmop(k)); % function to rootfind on
          if v, fprintf('rootfind in [%.16g, %.16g], N=%d:\n', k, ktop, N); end
          [kl el y u] = utils.intervalrootsboyd(f, [k ktop], io);
          [kl i] = sort(kl); el = el(i);
          if v, fprintf('found %d eigenvalues w/ %d fredholm det evals (ratio %.3g)\n', numel(kl), numel(y), numel(y)/numel(kl)); end
          kj = [kj; kl]; ej = [ej; el];
          k = ktop;                        % increment k and repeat
        end
        p.err.ej = ej;  % rootfinding error estimates (imaginary parts)
        p.kj = kj;
        
        if wantmodes
          if v, disp('computing eigenmodes at each eigenwavenumber...'); end
          p.solvemodescoeffs(meth, o); end  % compute all eigenfuncs
        
      else, error(fprintf('unknown method %s', meth));
      end
    
      if nargout>2, coj = p.coj; end        % outputs if needed...
      if nargout>3, ndj = p.ndj; end
    end
    
    function A = fillfredholmop(p, k) % .... helper for solvespectrum meth='fd'
    % FILLFREDHOLMOP - fill l^2 normed (sqrtwei) layer-potential matrix at k
      p.setoverallwavenumber(k);
      A = p.fillbcmatrix;   % fill matrix mapping density values to field values
      N = size(A,1); if size(A,2)~=N, error('A must be square!'); end
      % rescale so l^2 norms represent L^2 norms, as Jim Bremer also does...
      A = A .* repmat(1./p.sqrtwei, [N 1]); % note sqrtwei already left-applied 
      A = -2*A;     % for Dirichlet, -2 turns (D-1/2) into (1-2D) for det
    end
    
    function [coj ndj minsigj] = solvemodescoeffs(p, meth, o) %all evals in kwin
    % SOLVEMODESCOEFFS - compute modes and n-derivs at set of eigenwavenumbers
    %
    % [coj ndj minsigj] = SOLVEMODESCOEFFS(p, meth) takes eigenwavenumbers kj in
    %   p evp object and computes a mode at each. For degeneracies, eigenspaces
    %   are not computed. For meth='fd' (Fredholm det, the default),
    %   minsigj give the value of the minimum singular value of 1/2-D, which
    %   is a good error measure.
    %   This also verifies the list of kj as good approximate eigenwavenumbers.
    %   Data is saved in p object, but may also be passed out for convenience:
    %   Outputs: coj is stack of column vecs of coefficients of each mode,
    %   ndj - stack of col vecs giving boundary function data for each mode.
    %
    % [...] = SOLVEMODESCOEFFS(p, meth, opts) controls options such as:
    %   opts.verb = 0,1,... : verbosity (0=silent, 1=diagnostic text, ...)
    %
    % See also: EVP.SOLVESPECTRUM
      if nargin<3 | isempty(meth), meth='fd'; end   % default method
      if nargin<4, o = []; end
      if ~isfield(o, 'verb'), o.verb = 1; end; v = o.verb;  % default verbosity
      if numel(p.kj)<1, warning('no eigenwavenumbers kj in evp object!'); end
      p.coj = []; p.ndj = []; p.err.minsigj = [];  % erase any existing data
      for i=1:numel(p.kj), k = p.kj(i);
        [co nd e] = solvemodecoeffs(p, k, meth, o);
        p.coj = [p.coj co]; p.ndj = [p.ndj nd];  % stack (todo: realloc?)
        p.err.minsigj = [p.err.minsigj; e];
        if v, fprintf('\t mode #%d, k=%.16g: min sing val=%.3g\n',i,k,e); end
      end
      if nargout>=1, coj = p.coj; end        % outputs if needed...
      if nargout>=2, ndj = p.ndj; end
      if nargout>=3, minsigj = p.err.minsigj; end
    end
    
    function [co nd e] = solvemodecoeffs(p, k, meth, o) % ........ single efunc
    % SOLVEMODECOEFFS - helper (see SOLVEMODESCOEFFS) finds single efunc, O(N^3)
    %
    % see code for documentation
      if nargin<3 | isempty(meth), meth='fd'; end   % default method
      if nargin<4, o = []; end
      if ~isfield(o, 'iter'), o.iter = 1; end   % default iterative
      
      A = p.fillfredholmop(k);      % (I - 2D) in Dirichlet case
      if o.iter
        [u e v info] = utils.minsingvalvecs(A);  % iterative O(N^3) small const
        if info.flag, error('info.flag failed from utils.minsingvalvec!'); end
      else
        [U S V] = svd(A); u = U(:,end); v = V(:,end); e = S(end,end); % slow
      end
      co = v ./ p.sqrtwei.';        % density, convert from l^2 to value vector
      % L sv's are R sv's w/ D^*, so u gives boundary function...
      s = p.segs;                   % assume only a single segment!
      xdn = real(conj(s.x).*s.nx); w = sqrt(xdn/2); % sqrt Rellich bdry wei
      u = u/max(u);                 % make real-valued (imag part is error est)
      u = k * u/norm(w.*u);         % Rellich-normalize
      nd = u ./ p.sqrtwei.';        % convert from l^2 to value vector        
    end
    
    function [uj gx gy di js] = showmodes(p, o) % .............. plot all modes
    % SHOWMODES - compute and plot eigenmodes given their eigenvalues and coeffs
    %
    % SHOWMODES(p) computes and plots all eigenfunctions (modes) given by basis
    %   coefficients p.coj and eigenwavenumbers p.kj stored in evp problem
    %   object p. Modes are plotted in rectangular grid of subplots.
    %
    % [uj gx gy di js] = SHOWMODES(p) also outputs 3d array of real-valued mode
    %   values evaluated on a grid with x- and y-grids gx and gy, and where di
    %   is the array of inside-domain indices (as in problem.showsolution).
    %   js is the list of indices (ie wrt p.kj) that were computed in uj.
    %
    % SHOWMODES(p, opts) allows control of certain options including:
    %   opts.dx, opts.bb, opts.bdry - grid spacing, bounding-box, and whether
    %               to plot boundary, as in problem.showsolution.
    %   opts.eval - 'basis' uses existing basis for evaluation, 'grf' uses
    %               Green's rep. formula on p.ndj (Dirichlet BC only) (default)
    %   opts.kwin - only computes modes in given wavenumber window
    %   opts.inds - only computes modes in given (unordered) index set
    %
    % Notes/issues:
    %  * Normalization hasn't been considered much. For GRF case, they are
    %    correctly L2-normalized over the domain. Phase removal is a bad hack!
    %  * Degeneracies!
    %  * No account taken of non-simply connected domains (needs GRF for
    %    this to be set up...)
    %  * Neumann, Robin ... ?
    %  * Close-evaluation (J-expansion projection) for layer-potentials...?
    %
    % See also: EVP, PROBLEM.SHOWSOLUTION
      if nargin<2, o = []; end
      if ~isfield(o, 'eval'), o.eval='grf'; end, grf = strcmp(o.eval,'grf');
      if ~isfield(o, 'kwin'), o.kwin=p.kwin; end
      if ~isfield(o, 'inds'), o.inds=1:numel(p.kj); end
      if ~isfield(o, 'bdry'), o.bdry = 0; end
      wantdata = nargout>0;
      inlist = 0*p.kj; inlist(o.inds) = 1;         % true if in index list
      js = find(p.kj>o.kwin(1) & p.kj<o.kwin(2) & inlist); % indices to plot
      if isempty(js), warning('no modes found in evp, or list or window!'); end
      
      n = numel(js);
      figure; nac = ceil(sqrt(n)); ndn = ceil(n/nac); % # subplots across & down
      if grf, d = utils.copy(p.doms);  % set up for GRF evaluations in loop...
        d.clearbases; d.addlayerpot(d.seg, 's');
        pe = bvp(d); pe.setupbasisdofs; % temporary new problem instance
      else pe = p; end                  % just use existing problem
        
      for i=1:n, j = js(i);      % ---- loop over selected eigenvalues
        pe.setoverallwavenumber(p.kj(j));                         % get this k
        if grf, pe.co = p.ndj(:,j); else, pe.co = p.coj(:,j); end % get coeffs
        [u gx gy di] = pe.gridsolution(o); % either conventional or GRF evalu
        if ~grf  % rotate phase of u (est at int pt) so that essentially real
          uint = u(floor(numel(gx)/2),floor(numel(gy)/2)); % HACK guess int pt!
          u = u * conj(uint/abs(uint));  % unphase u
        end
        u = real(u);                     % saves half the memory (real-valued)
        if wantdata
          if i==1, uj = nan(size(u,1), size(u,2), n); end % allocate output
          uj(:,:,i) = u;                                  % copy data to array
        end
        tsubplot(ndn, nac, i); imagesc(gx, gy, u, 'alphadata', ~isnan(di));
        caxis(3.5*[-1 1]/sqrt(p.doms.area)); % rescale values based on area
        axis equal tight; colormap(jet(256)); axis off;
        set(gca,'ydir','normal'); hold on; if o.bdry, pr.showbdry; end
        drawnow;                                          % for fun
      end                        % ----
    end
    
    function weylsmoothcheck(p, o)
      % ...
      % this will dump kj.^2 as Gaussians into regular E array, compare to Weyl
    end
  end
  
  % --------------------------------------------------------------------
  methods(Static)    % these don't need EVP obj to exist to call them...
    [k_weyl] = weylcountcheck(k_lo, k, perim, area, dkfig)   % old Weyl routine
  end
end
