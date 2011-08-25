% EVP - define a Laplace eigenvalue problem for a single interior domain
%
% p = evp(dom) creates an EVP object from the domain dom, comprising the
%   domain along with its segments the boundary conditions they carry, and
%   the basis sets associated with the domain.
%
%   Preliminary version: only one domain is supported, with Dirichlet or
%   Neumann BCs. However, multiple segments and multiply-connected seem to
%   work ok.
%
%   Solution methods available are:
%    * root-search on the Fredholm determinant of Id -+ 2D (double-layer op).
%    * min-singular-value root search on (Id -+ 2D).
%    * crude version of MPS on a single domain with multiple segments.
%    * weighted-NtD scaling method with various prediction orders

% Copyright (C) 2010 - 2011, Alex Barnett, Timo Betcke

classdef evp < problem & handle
  properties
    kj                        % eigenwavenumbers = sqrt(eigenvalues) (col vec)
    coj                       % basis coeffs (columns of array) of modes
    err                       % struct for error estimates errors of kj:
                              %   err.ej - rootfinding errors (col vec)
                              %   err.minsigj - min op sing vals (col vec)
    ndj                       % normal derivatives (cols of array) of modes
    kwin                      % wavenumber window in which kj were requested
    intpts                    % interior pointset (used only for MPS)
  end

  methods % ----------------------------------------- particular to EVP
    function pr = evp(doms) % .........................constructor
      if nargin==0, return; end    % needs empty constructor
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
    % [kj] = SOLVESPECTRUM(p, [klo khi]) uses default 'fd' method to find all
    %   eigenwavenumbers kj lying in the wavenumber interval [klo, khi], for the
    %   eigenvalue problem object p.
    %   Laplacian eigenvalues are simply the squares of the kj.
    %   The problem also stores the kj list as p.kj
    %
    % [kj err coj ndj] = ... also returns basis coeffcients, normal derivatives,
    %   and error estimates struct, if possible. Certain methods may return
    %   empty arrays for some of these.
    %
    % SOLVESPECTRUM(p, [klo khi]) with no output argument stores kj and err
    %   (and if opts.modes=1, also ndj) as property fields in the object p.
    %
    % [...] = SOLVESPECTRUM(p, [klo khi], meth, opts) allows one to choose the
    %   method,
    %     meth = 'fd' : Fredholm det root search (layer-potential bases only)
    %                   Error estimates are dist to real axis of Boyd roots.
    %                   Id +- 2(double-layer operator), for a single domain, w/
    %                   Dirichlet or Neumann BCs everywhere on the boundary.
    %     meth = 'ntd' : weighted-Neumann-to-Dirichlet scaling method
    %                   (star-shaped domains with single segment only)
    %                   coj is returned empty since all bdry data is in ndj.
    %     meth = 'ms' : iterative search for minimum singular value of Id+-2D.
    %                   (considered state-of-art in literature).
    %
    %   and/or choose options,
    %     opts.verb = 0,1,... : verbosity (0=silent, 1=diagnostic text, ...)
    %     opts.modes = 0,1 : if true, compute modes (storing in p)
    %                  (note, for ntd khat='o', requesting modes not extra work)
    %     opts.eps : size of k-window over which meth='ntd' extrapolates
    %                (default is 0.2/(max radius of domain))
    %     opts.khat = 'ntd' k_hat predict method: 'l' linear, 'o' frozen-f ODE.
    %     opts.fhat = 'ntd' f_hat predict method: 'f' use f*, 'l' lin, 's' 2nd-o
    %
    % Notes / issues:
    % * 'fd' and 'ms' require p to have appropriate layerpot basis on domain
    %   ('ntd' doesn't need this)
    % * Still very alpha code! Single domain, but can have multiple segs.
    % * keep more data of all the det evals, pass out?
      if nargin<3 | isempty(meth), meth='fd'; end   % default method
      if nargin<4, o = []; end
      if ~isfield(o, 'tol'), o.tol = 1e-8; end   % tolerance for all rootfinds
      if ~isfield(o, 'verb'), o.verb = 1; end; v = o.verb; % default verbosity
      if ~isfield(o, 'modes'), o.modes = 0; end
      if numel(p.doms)~=1, warning('problem has more than one domain...'); end
      d = p.doms;                            % the domain(s)
      d1 = d(1);            % the first, and we hope only, domain
      wantmodes = nargout>2 | o.modes;
      if v, fprintf('solvespectrum: wantmodes = %d\n', wantmodes); end
      klo = kwin(1); khi = kwin(2); p.kwin = kwin;
      p.fillquadwei; M = numel(p.sqrtwei);   % total # boundary pts
      if v, fprintf('mean M ppw @ khi = %.3g\n', 2*pi*M/khi/d1.perim); end
      
      if strcmp(meth,'fd') % ............... Fredholm det method..............
        k = klo; kj = []; ej = []; N = M;  % square system
        kscale = 70/d1.diam;             % est k at which level density ~ 20
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
        
      elseif strcmp(meth,'ntd') % ............... NtD scaling method...........
        if ~isfield(o, 'eps'), o.eps = 0.2/d1.diam; end % defaults
        if ~isfield(o, 'khat'), o.khat = 'o'; end         % (most acc)
        if ~isfield(o, 'fhat'), o.fhat = 's'; end         % (")
        sx = d1.x; snx = d1.nx; sp = d1.speed; % quad info for the one domain
        xn = real(conj(sx).*snx); ww = p.sqrtwei.'.^2./xn; % x.n bdry func
        ixnip = @(f,g) sum(ww.*conj(f).*g);    % 1/(x.n)-weighted inner prod
        ixnrm = @(g) sqrt(real(ixnip(g,g)));   % 1/(x.n)-weighted bdry norm
        xt = real(conj(sx).*1i.*snx);
        if wantmodes | ~strcmp(o.khat,'l')
          D = 2*pi*repmat(1./sp,[1 M]) .* circulant(quadr.perispecdiffrow(M));
          xnt = D*xn; xtt = D*xt; m = (xn.*xtt-xt.*xnt)./xn;  % Andrew's m func
        end  % note avoided storing and computing V matrix
        k = klo; kj = []; err.imbej = []; err.koj = []; p.ndj = []; % prealloc?
        while k<khi      % ==== loop over windows
          ktop = min(k + o.eps, khi); % top of current k window
          if v, fprintf('NtD k* = %.14f ...\t', k); end
          if ~wantmodes & strcmp(o.khat,'l') % eigenvectors not needed
            d = p.NtDspectrum(k, o);
          else, [d U] = p.NtDspectrum(k, o); end
          betamin = -1.1*o.eps/k;  % little penalty for excess factor here
          % loop over all negative eigvals and keep those khats in window...
          rd = real(d)'; % note rd is row vec for for-loop find work
          kl = []; imbel = []; ndl = []; kol = [];
          for i=find(rd>betamin & rd<0), be = rd(i);
            if exist('U','var'), f = U(:,i); end  % NtD eigenfunc at kstar
            khat = k/(1+be);  % basic eigenfreq approximation
            if strcmp(o.khat, 'o')
              Df = D*f; mf = m.*f; % compute ODE constants...
              A = ixnrm(xn.*Df)^2; B = real(ixnip(f,mf)); C = ixnrm(xn.*f)^2;
              khat = evp.solvebetaode(k, be, A, B, C); % frozen A,B,C
            end
            if wantmodes                     % recon boundary efunc...
              if strcmp(o.fhat, 'f'), fhat = f;
              else         % use 1st or 2nd deriv of f
                if ~exist('Df','var'), Df = D*f; mf = m.*f; end 
                Vf = xt.*Df; Vff = ixnip(f,Vf); dfdk = (Vf + mf + Vff*f)/k;
                fhat = f + (khat-k)*dfdk; % 1st, uses best avail khat. Normed
                if strcmp(o.fhat, 's')       % use 2nd-deriv of f:
                % Alex formula, leading order bits only, no c parallel cmpt...
                  ddfdkk = xn.*xn.*f + (xt.*(D*(Vf+mf)) + m.*(Vf+mf) + ...
                                        xn.*(D*(xn.*(Df))))/k^2;
                  fhat = fhat + (khat-k)^2*ddfdkk/2;
                end
              end
              fhat = fhat / ixnrm(fhat);   % normalize (avoids needing c)
              nd = sqrt(2) * khat * fhat ./ xn; % normal-deriv function 
            else nd = []; end
            if khat>=k & khat<ktop   % if in window, keep it
              kl = [kl; khat]; imbel = [imbel; imag(d(i))]; % keep Im(beta)
              ndl = [ndl nd]; kol = [kol; k];
            end
          end
          if v, fprintf('found %d eigvals in [%g,%g]\n', numel(kl),k,ktop); end
          [dummy i] = sort(kl); % reorder all the l-arrays when append to lists
          kj = [kj; kl(i)]; err.imbej = [err.imbej; imbel(i)];
          err.koj = [err.koj; kol(i)];
          if wantmodes, p.ndj = [p.ndj ndl(:,i)]; end
          k = ktop;                        % increment k and repeat
        end            % ==== end window loop
        p.kj = kj; p.err = err; p.coj = [];
        
      elseif strcmp(meth,'ms') % ........... SVD iter min search method........
        io = o; io.xtol = o.tol;
        dE = 4*pi/d1.area;  % Weyl mean level spacing in E=k^2
        ppls = 5;          % user param: # grid points per mean levelspacing
        ng = ceil(ppls * (khi^2-klo^2)/dE + 1);
        if v, fprintf('# k-gridpts = %d. Doing gridminfit...\n', ng); end
        g = sqrt(linspace(klo^2, khi^2, ng));
        f = @(k) svd(p.fillfredholmop(k)); % vector function to rootfind on
        [kj sj p.err.mininfo] = evp.gridminfit(f, g, io); % the meat: do search
        % currently discards higher sing vals - keep them?
        p.err.ej = sqrt(min(sj,[],1));  % min sing val error estimates
        kj = kj(:); p.err.ej = p.err.ej(:); % make col vecs
        p.kj = kj;
        
        if wantmodes  % same as in 'fd'
          if v, disp('computing eigenmodes at each eigenwavenumber...'); end
          p.solvemodescoeffs(meth, o); end  % compute all eigenfuncs

      else, error(fprintf('unknown method %s', meth)); % ...................
      end
    
      if nargout>2, coj = p.coj; end    % outputs if needed (might be empty)...
      if nargout>3, ndj = p.ndj; end
    end
    
    function [t co] = tension(p, k, o) % crude, one method for now
      if nargin<3, o = []; end
      if ~isfield(o, 'eps'), o.eps = 1e-14; end
      wantvec = nargout>1;
      p.setoverallwavenumber(k);
      A = p.fillbcmatrix;
      B = p.evalbases(p.intpts);
      if 0 %    unregularized GSVD
        if wantvec
          [UU VV X C S] = gsvd(A,B,0);  % co gives eigvecs as rows
          co = X;
          t = sqrt(diag(C'*C)./diag(S'*S));
        else, t = gsvd(A,B); end
      else
      % reg GSVD. Timo's code from hassell/bnds/tensionsq.m :
      [Q,R]=qr([A;B],0); [U,S,V]=svd(R); S=diag(S);
      ii = abs(S)>o.eps*max(abs(S));   % note I scaled it to max(sig val)
      %numel(ii) % rank
      Q=Q*U(:,ii);  % cols of Q now onb for Col[A;B] at numerical rank.
      if wantvec
        [UU VV X C S] = gsvd(Q(1:size(A,1),:),Q(size(A,1)+1:end,:));
        co = R*(V(:,ii) * X);   % rotate eigvecs back: coeffs wrong ???
        t = sqrt(diag(C'*C)./diag(S'*S));
      else
        t = gsvd(Q(1:size(A,1),:),Q(size(A,1)+1:end,:));
      end
      end
    end
    
    function [kj err coj sweep] = crudempssolvespectrum(p, kwin) % for Nilima
      dk = (kwin(2)-kwin(1))/100;
      ks = kwin(1):dk:kwin(2); ts = nan(p.N, numel(ks));   % tension sweep
      disp('sweeping...')
      for i=1:numel(ks), k=ks(i);
        t = p.tension(k); ts(1:numel(t),i) = t;
      end
      sweep.ts = ts; sweep.ks = ks;
               figure; plot(ks, ts, '-+'); axis([kwin(1) kwin(2) 0 1]); drawnow; kj =0; err=0; coj = 0; return;
      mt = min(ts, [], 1);
      kinit = ks(find(mt<0.01));  % initial guesses
      op = optimset('tolX',1e-12,'display', 'iter');
      disp('optimizing...')
      for j=1:numel(kinit), kinit(j)
        [kj(j) tj(j) flag(j) outp(j)] = fminbnd(@(k) min(p.tension(k))^2, ...
                                                kinit(j)-dk,kinit(j)+dk,op);
      end
      [kj,I] = sort(kj,'ascend');
      newj = [1 1+find(abs(diff(kj))>1e-7)]; % remove duplicates: need max k err
      kj = kj(newj);
      coj = nan(p.N, numel(kj));  % now get coeffs... (I don't believe these)
      disp('getting coeffs...')
      for j=1:numel(kj)
        [t co] = p.tension(kj(j)); coj(:,j) = co(:,1); err.tj(j) = t(1);
      end
    end
    
    
    function setupintpts(p, I) % crappy slow for now. 3/24/11
      intpts = zeros(I,1);
      go = p.gridboundingbox; bb = go.bb; clear go
      i = 1;
      while i<I
        intpts(i) = bb(1)+rand(1)*(bb(2)-bb(1)) + ...
            1i*(bb(3)+rand(1)*(bb(4)-bb(3)));
        if ~isnan(p.domainindices(pointset(intpts(i)))), i=i+1; end
      end
      p.intpts = pointset(intpts);
    end
    
    function A = fillfredholmop(p, k) % .... helper for solvespectrum meth='fd'
    % FILLFREDHOLMOP - fill l^2 normed (sqrtwei) layer-potential matrix at k
    %
    % This works for Dirichlet (bc='D', lp='d') or Neumann (bc='N', lp='s')
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
      if nargin<2 | isempty(meth), meth='fd'; end   % default method
      if nargin<3, o = []; end
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
    % see code for documentation. Dir is normalized correctly; Neu is not.
    % Braxton needed multiple segs, so done.
      if nargin<3 | isempty(meth), meth='fd'; end   % default method
      if nargin<4, o = []; end
      if ~isfield(o, 'iter'), o.iter = 1; end   % default iterative
      
      A = p.fillfredholmop(k);      % (I - 2D) in Dirichlet, (I + 2D^T) Neu
      if o.iter
        [u e v info] = utils.minsingvalvecs(A);  % iterative O(N^3) small const
        if info.flag, error('info.flag failed from utils.minsingvalvec; try setting opts.iter=0 instead!'); end
      else
        [U S V] = svd(A); u = U(:,end); v = V(:,end); e = S(end,end); % slow
      end
      if p.segs(1).a==0, z = u; u = v; v = z; end % for Neu case, swap u,v HACK
      % need to correctly normalize Neu case now... (specrowdiff needed)
      % ...
      co = v ./ p.sqrtwei.';        % density, convert from l^2 to value vector
      % L sv's are R sv's w/ D^*, so u gives boundary function...
      x = vertcat(p.segs.x); nx = vertcat(p.segs.nx); % allow multiple segs
      xdn = real(conj(x).*nx); w = sqrt(xdn/2); % sqrt Rellich bdry wei
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
      if isempty(js)
        warning('no eigenvalues found in evp, or list or window!'); end      
      
      if grf
        if isempty(p.ndj), error('GRF method: no ndj bdry func info found!');end
        d = utils.copy(p.doms);  % set up for GRF evaluations in loop...
        d.clearbases; d.addlayerpot(d.seg, 's');
        pe = bvp(d); pe.setupbasisdofs; % temporary new problem instance
      else
        if isempty(p.coj), error('basis method: no coj coeffs info found!');end
        pe = p;                  % just use existing problem
      end
      n = numel(js);
      figure; nac = ceil(sqrt(n)); ndn = ceil(n/nac); % # subplots across & down

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
    
    % methods needing an EVP object, but defined by separate files...
    [d V] = NtDspectrum(p, kstar, opts)
   end % methods
  
  % --------------------------------------------------------------------
  methods(Static)    % these don't need EVP obj to exist to call them...
    [k_weyl] = weylcountcheck(k_lo, k, perim, area, dkfig)   % old Weyl routine
    kh = solvebetaode(kstar, beo, A, B, C, Ap, Bp, Cp)  % NtD scaling helper
    [xm ym info] = gridminfit(f, g, o)   % gridded vectorial minimizer
    [xm fm info] = iterparabolafit(f, x, y, opts)   % minimization helper
    [A,B,C] = para_fit(e, f)        % parabolic helper
  end
end
