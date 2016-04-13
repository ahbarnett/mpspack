% EVP - define a Laplace eigenvalue problem for a single interior domain
%
% p = evp(dom) creates an EVP object from the domain dom, comprising the
%   domain along with its segments the boundary conditions they carry, and
%   the basis sets associated with the domain.
%
%   Preliminary version: only one domain is supported, with Dirichlet or
%   Neumann BCs. However, multiple segments and multiply-connected seem to
%   work ok, assuming the method is compatible with corners.
%
%   Solution methods available are:
%    * root-search on the Fredholm determinant of Id -+ 2D (double-layer op).
%    * min-singular-value root search on (Id -+ 2D).
%    * crude version of MPS on a single domain with multiple segments.
%    * weighted-NtD scaling method (Dir BCs) with various prediction orders
%    * filtered-DtN scaling method (Neu BCs) w/ various prediction orders
%      (for Neumann MPS inclusion methods of 2015 see examples/neumann_evp)

% Copyright (C) 2010 - 2016, Alex Barnett

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
  
    function [kj err coj ndj dat] = solvespectrum(p, kwin, meth, o)
    % SOLVESPECTRUM - compute spectrum, and possibly modes, in wavenumber window
    %
    % [kj] = SOLVESPECTRUM(p, [klo khi]) uses default 'fd' method to find all
    %   eigenwavenumbers kj lying in the wavenumber interval [klo, khi], for the
    %   eigenvalue problem object p. Everywhere Dirichlet or Neumann supported
    %   in all methods apart fron 'ntd'.
    %   Laplacian eigenvalues are simply the squares of the kj.
    %   The problem also stores the kj list as p.kj
    %
    % [kj err coj ndj] = ... also returns basis coeffcients (i.e. the relevant
    %   layer density function: tau for Dir BCs, sigma for Neu), the boundary
    %   data (normal derivatives for Dir, values for Neu),
    %   and an error estimate struct, if possible. Certain methods may return
    %   empty arrays for some of these.
    %
    % [kj err coj ndj dat] = ... also returns NtD spectral sweep data object
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
    %     meth = 'int' : interpolating wNtD or spectrally-weighted DtN
    %                   (star-shaped domains with single segment only)
    %                   coj is returned empty since all bdry data is in ndj.
    %     meth = 'ms' : iterative search for minimum singular value of Id+-2D.
    %                   (considered state-of-art in literature).
    %
    %   and/or choose options,
    %     opts.verb = 0,1,... : verbosity (0=silent, 1=diagnostic text, ...)
    %     opts.modes = 0,1 : if true, compute modes (storing in p)
    %                  (note, for ntd khat='o', requesting modes no extra work)
    %     opts.eps : size of k-window over which meth='ntd' extrapolates
    %                (default is 0.2/(max radius of domain))
    %     opts.khat = 'ntd' k_hat predict method: 'l' linear, 'o' frozen-f ODE.
    %                 'n' Neumann linear, 'p' Neumann frozen ODE.
    %     opts.p    = 2,4,.. Lagrange interpolation order for meth='int'.
    %     opts.itype = 'l' Lagrange, 'h' Hermite, for beta(k) and f(k), etc.
    %     opts.fhat = 'ntd' f_hat predict method: 'f' use f*, 'l' lin, 's' 2nd-o
    %     opts.dat : if nonempty, reuse as 'ntd' or 'int' spectral sweep data.
    %
    % Notes / issues:
    % * 'fd' and 'ms' require p to have appropriate layerpot basis on domain
    %   and BCS on segment, already set up ('ntd' doesn't need this)
    % * Still very alpha code! Single domain, but can have multiple segs.
    % * keep more data of all the det evals, pass out?
      if nargin<3 | isempty(meth), meth='fd'; end   % default method
      if nargin<4, o = []; end
      if ~isfield(o, 'tol'), o.tol = 1e-8; end   % tolerance for all rootfinds
      if ~isfield(o, 'verb'), o.verb = 1; end; v = o.verb; % default verbosity
      wantmodes = nargout>2; if isfield(o, 'modes'), wantmodes = o.modes; end
      datavail=0; dat = []; if isfield(o, 'dat'), dat = o.dat; datavail=1; end
      outdat = nargout>4;       % want saving of spectral data?
      if numel(p.doms)~=1, warning('problem has more than one domain...'); end
      d = p.doms;                            % the domain(s)
      d1 = d(1);            % the first, and for some methds only, domain
      sx = d1.x; snx = d1.nx; sp = d1.speed; % quad info for the one domain
      ip = @(f,g) sum(p.sqrtwei(:).^2.*conj(f(:)).*g(:));  % unwei inner prod
      xn = real(conj(sx).*snx);                  % x.n
      if strcmp(meth,'int') || strcmp(meth,'ntd')
        ww = p.sqrtwei.'.^2./xn; % 1/x.n bdry func w/ weights
        ixnip = @(f,g) sum(ww(:).*conj(f(:)).*g(:)); % 1/(x.n)-wei inner prod
        ixnrm = @(g) sqrt(real(ixnip(g,g)));   % 1/(x.n)-weighted bdry norm
      end
      if v, fprintf('solvespectrum: wantmodes = %d\n', wantmodes); end
      klo = kwin(1); khi = kwin(2); p.kwin = kwin;
      p.fillquadwei; M = numel(p.sqrtwei);   % total # boundary pts
      if v, fprintf('mean M ppw @ khi = %.3g\n', 2*pi*M/khi/d1.perim); end
      neu = (d1.seg(1).a==0);                % is 1 seg Neumann BCs?
      
      if strcmp(meth,'fd') % ....... Fredholm det Boyd rootfind method..........
        k = klo; kj = []; ej = []; N = M;  % square system
        kscale = 70/d1.diam;             % est k at which level density ~ 20
        io.Ftol = 1e-10; io.tol = o.tol;   % allows leeway in collecting roots
        io.disp = v;                       % opts for interval search
        y = []; u = [];         % keep fred det eval data
        while k<khi
          ktop = min(k + kscale*min(1/k,0.2), khi); % top of current k window
          % option to reset N here, requadrature all segs (in proportion?)
          % ...
          %f = @(k) det(eye(N) - 2*layerpot.D(k, s)); % for one closed segment
          f = @(k) det(p.fillfredholmop(k)); % function to rootfind on
          if v, fprintf('rootfind in [%.16g, %.16g], N=%d:\n', k, ktop, N); end
          [kl el yl ul] = utils.intervalrootsboyd(f, [k ktop], io);
          [kl i] = sort(kl); el = el(i);
          if v, fprintf('found %d eigenvalues w/ %d fredholm det evals (ratio %.3g)\n', numel(kl), numel(y), numel(y)/numel(kl)); end
          kj = [kj; kl]; ej = [ej; el];
          y = [y yl(end:-1:1)]; u = [u ul(end:-1:1)]; % make in incr k order
          k = ktop;                        % increment k and repeat
        end
        p.err.ej = ej;  % rootfinding error estimates (imaginary parts)
        p.err.y = y; p.err.u = u;  % output fd data
        p.kj = kj;
        
        if wantmodes
          if v, disp('computing eigenmodes at each eigenwavenumber...'); end
          p.solvemodescoeffs('ms', o); end  % compute all eigenfuncs via ms
        
      elseif strcmp(meth,'int') % ......... interpolate wNtD or wDtN flow.......
        % Note the f-inner-prod matching will fail for multiple EVs - ** fix.
        if ~isfield(o, 'eps'), o.eps = 0.2/d1.diam; end % defaults
        if ~isfield(o, 'p'), o.p = 2; end % default interpolation order
        if ~isfield(o,'itype'), o.itype = 'l'; end
        nspill = ceil(o.p/2)-1; % how many extra k samples over the ends
        clear filteredDtNspectrum  % since persistent filter mats used in it
        nw = ceil((khi-klo)/o.eps); o.eps = (khi-klo)/nw;       % # windows
        if ~datavail  % fill whole spectral dat... (the expensive part)
          dat.nos = -nspill:nw+nspill;    % indices of k slices spectral data
          if v, fprintf('meth=int make dat: nw=%d nspill=%d eps=%g\n',nw,nspill,o.eps); end
          dat.kos = klo + o.eps * dat.nos;      % k's where spectra computed
          for n=1:numel(dat.nos), k = dat.kos(n); % ----- dat loop over windows
            if v, fprintf('wDtN int dat k* = %.15g ...    \t', k); end
            [d U] = p.filteredDtNspectrum(k, o);
            betakeep = 1.3*o.eps*k;  % since beta_dot is roughly -k
            ii = find(abs(real(d))<=betakeep*o.p/2);      % indices to keep
            [~,j] = sort(real(d(ii)));        % beta's in incr real part
            if v, fprintf('%d EVs of wDtN kept\n', numel(ii)); end
            bes = d(ii(j)); fs = U(:,ii(j)); % beta's (assume ascending), f's
            dat.bes(n,1:numel(bes)) = bes;   % save to dat
            dat.nes(n) = numel(bes);         % number of eigenvalues kept
            for i=1:numel(bes), dat.fs(n,:,i) = fs(:,i)/ixnrm(fs(:,i)); end
          end                                     % -------
        end
        kj = []; err.imbej = []; err.koj = []; p.ndj = []; p.coj = []; % koj will store nearest slice to each kpred
        if datavail % now use dat to do the prediction...
        i = find(dat.nos==0);          % moving index into dat objects
        for n=0:nw-1                             % ----- pred loop over windows
          k0 = dat.kos(i);           % bottom end of k window
          %plot(k0, real(dat.bes(i,:))./(dat.bes(i,:)~=0), 'k.'); hold on; %dat.bes(i,:), dat.bes(i+1,:)
          betakeep = 1.1*o.eps*k0;   % smaller fudge factor in central window
          jj=find(dat.bes(i,:)>0 & dat.bes(i,:)<betakeep);   % indices of beta's which flow could give kpred in this window [k0,k0+eps] (beta=0 is bogus)
          if numel(jj)>0     % if no beta>0 then nothing to predict in window
            q(1+nspill) = jj(1);     % index in current slice of smallest beta>0
            fj = dat.fs(i,:,q(1+nspill));
            for l=[-nspill:-1, 1:nspill+1];  % loop over lower & upper k slices
              %plot(dat.kos(i+l), real(dat.bes(i+l,:))./(dat.bes(i+l,:)~=0), 'b.');
              m = dat.nes(i+l);   % how many beta's to check at this slice
              ips = nan(m,1);   % inner-products
              for j=1:m, ips(j) = ixnip(fj,dat.fs(i+l,:,j)); end
              jyes = find(abs(ips)>0.9);   % indices which match efunc well
              if numel(jyes)~=1, warning(sprintf('j match prob! k=%g l=%d numel(jyes)=%d',k0,l,numel(jyes))); ips, end
              if ips(jyes)<0, dat.fs(i+l,:,jyes) = -dat.fs(i+l,:,jyes); end
              q(1+nspill+l) = jyes;  % which associates w/ jj(1) in this slice
            end  % now each index q(1:2*(nspill+1)) should have been assigned
            for j=jj
              x = -nspill:1+nspill; % poly fit ordinates (in local window-lens)
              y = []; for l=1:2*(nspill+1)      % poly fit data
                y(l) = real(dat.bes(i+l-1-nspill,q(l))); end
              pp = polyfit(x,y,o.p-1); x0 = fzero(@(x) polyval(pp,x),0.5);
              kpred = k0 + x0*o.eps;
              if kpred>k0 && kpred<=k0+o.eps  % collect the eigenpair
                plot(k0+o.eps*x,y,'-'); plot(real(kpred),0,'r+');
                
                t = x0;   % fix up f interp to use lagrange??? *******
                fpred = t*dat.fs(i+1,:,q(1+nspill)) + (1-t)*dat.fs(i,:,q(2+nspill)); % NB f's signs were fixed earlier so interp makes sense
                
                kj = [kj;kpred]; p.ndj = [p.ndj fpred]; % *** convert f to u!
                if t<0.5; err.koj = [err.koj;k0];  % save nearest kstar as koj
                else err.koj = [err.koj;k0+o.eps]; end
              end
              q = q + 1;   % increment all indices in fitting slices together
              if j~=jj(end) && dat.bes(i+1,q(2+nspill))>0, break, end  % not in window
            end            
            %          for l=2+nspill:o.p     % find starting indices q at larger k slices
          end
          i = i+1;
        end                                      % ------
        p.kj = kj; p.err = err;
        end % pred step
        
        
      elseif strcmp(meth,'ntd') % .......... NtD "scaling" type method..........
        if ~isfield(o, 'eps'), o.eps = 0.2/d1.diam; end % defaults
        if ~isfield(o, 'khat'), if neu, o.khat = 'p'; else, o.khat = 'r'; end, end        % (most acc)
        if ~isfield(o, 'fhat'), o.fhat = 's'; end         % (")
        xt = real(conj(sx).*1i.*snx); % build geom funcs upon xn, etc
        if wantmodes | ~strcmp(o.khat,'l') % todo: replace D by FFT application
          D = 2*pi*repmat(1./sp,[1 M]) .* circulant(quadr.perispecdiffrow(M));
          xnt = D*xn; xtt = D*xt; m = (xn.*xtt-xt.*xnt)./xn;  % Andrew's m func
        end  % note avoided storing and computing V matrix
        kbot = klo; % starting value for bottom of interval - what's looped over
        clear filteredDtNspectrum  % since persistent filter mats used in it
        kj = []; err.imbej = []; err.koj = []; p.ndj = []; % prealloc?
        n = 1;    % index to write to or read out of, in dat obj
        while kbot<khi      % ==== loop over windows (kstart is called just k)
          ktop = min(kbot + o.eps, khi); % top of current k window
          if ~datavail     % compute NtD spec and possibly eigenvectors
            if ~neu
              k = kbot; if v, fprintf('NtD k* = %.15g ...    \t', k); end
              if ~wantmodes & strcmp(o.khat,'l') % Dir; evecs not needed
              clear U; d = p.NtDspectrum(k, o);
                    else, [d U] = p.NtDspectrum(k, o); end
              betakeep = 1.1*o.eps/k;  % beta range to keep in dat (RAM hog)
            else, k = ktop;  % predict from top of interval
              if v, fprintf('DtN k* = %.15g ...\t', k); end
              if ~wantmodes % Neu; eigenvectors not needed
                clear U; d = p.filteredDtNspectrum(k, o); % use k=ktop
              else, [d U] = p.filteredDtNspectrum(k, o); end
              betakeep = 1.1*o.eps*k;  % since beta_dot is roughly -k
            end
            % loop over all negative eigvals and keep those khats in window...
            rd = real(d)';
            ii = find(rd>-betakeep & rd<=0);              % std eigvals to keep
            bes = d(ii); if exist('U','var'), fs = U(:,ii); end   % keep
            if outdat, dat.bes(n,1:numel(bes)) = bes;   % save to dat object
              dat.nes(n) = numel(bes);       % number of eigenvalues kept
              if exist('fs','var'), dat.fs(n,:,1:numel(bes)) = fs; end
            end
          else    % instead get bes, fs, from n'th spectral data object
            ne = dat.nes(n); bes = dat.bes(n,1:ne);
            if ~neu, k = kbot; else, k = ktop; end % recover the k* used
            if isfield(dat,'fs'), fs = reshape(dat.fs(n,:,1:ne),M,ne); end
          end
          kl = []; imbel = []; ndl = []; kol = [];
          for i=1:numel(bes), be = real(bes(i));  % only use real part for comp
            if exist('fs','var'), f = fs(:,i); f = f/ixnrm(f); end % f at kstar
            if ~neu          % Dir case
              khat = k/(1+be);  % basic linear eigenfreq approximation
            else
              khat = k + 1.000*be/k;  % basic linear eigenfreq approx
  %            khat = (k + sqrt(k^2+4*be))/2;  % basic linear eigenfreq approx
            end
            if strcmp(o.khat, 'r')   % Ricatti exact (Andrew)
              Df = D*f; mf = m.*f; % compute ODE constants...
              kavg = (khat+k)/2;  % kavg=k; % average k to use for Ricatti
              A = kavg^2*ixnrm(xn.*f)^2-ixnrm(xn.*Df)^2; B = -real(ixnip(f,mf));
              mu = sqrt(A - B^2/4); B2mu = B/(2*mu);
              khat = k*exp((atan(B2mu)-atan(A*be/mu + B2mu))/mu);
            elseif strcmp(o.khat, 'o')  % ODE numerical (slowish)
              Df = D*f; mf = m.*f; % compute ODE constants...
              %A = k^2*ixnrm(xn.*f)^2-ixnrm(xn.*Df)^2; B = -real(ixnip(f,mf));
              %khat = evp.solvebetaode(k, be, A, -B); % frozen A,B
              A = ixnrm(xn.*Df)^2; B = real(ixnip(f,mf)); C = ixnrm(xn.*f)^2;
              khat = evp.solvebetaode(k, be, A, B, C); % frozen A,B,C
              
            elseif strcmp(o.khat, 'n')  % Neumann slope expt fix 1/13/14
              ns = -M/4:M/4;  % code segment for iF from filteredDtNspectrum:
              L = sum(d1.w);
              arcl = real(ifft(-1i*fft(sp).*[0 1./(1:M/2-1) 0 -1./(M/2-1:-1:1)]')); % spectral approx to sampled cumulative arclength
              arcl = arcl'/2/pi + L*d1.seg(1).t'; % add growing component
              P = exp((-2i*pi/L)*ns'*arcl); % ns=freqs. Dense DFT mat
              Pt = P'; PW = P.*repmat(d1.w/L,[numel(ns) 1]);
              kfilt = k; xin = 2*pi/L/kfilt*ns;  % freq xi grid
              Fk = max(real(sqrt(1-xin.^2)),0.1*kfilt^(-1/3)); % match prefac !!
              iF = Pt*(diag(Fk.^-1)*PW); % note no high freq cutoff
              u = iF * (f./xn); % get neu bdry data from efunc
%              figure; plot(arcl, real([f u]), '+-'); % check iF does sensible
              dbe = sum(d1.w' .* xn .* u.^2);  % form for beta-dot...
              dbe = dbe - sum(d1.w' .* xn .* (D*u).^2)/k^2;
              dbe = -dbe*k; %dbe = -k*(1 + 0.03*rand(1)); % jiggle them
              khat = k - be/dbe;  % use as estimator
            
            elseif strcmp(o.khat, 'p')  % Neu ODE numerical expt 1/14/14
              ns = -M/4:M/4;  % code segment for iF from filteredDtNspectrum:
              L = sum(d1.w);
              arcl = real(ifft(-1i*fft(sp).*[0 1./(1:M/2-1) 0 -1./(M/2-1:-1:1)]')); % spectral approx to sampled cumulative arclength
              arcl = arcl'/2/pi + L*d1.seg(1).t'; % add growing component
              P = exp((-2i*pi/L)*ns'*arcl); % ns=freqs. Dense DFT mat
              Pt = P'; PW = P.*repmat(d1.w/L,[numel(ns) 1]);
              kfilt = k; xin = 2*pi/L/kfilt*ns;  % freq xi grid
              Fk = max(real(sqrt(1-xin.^2)),0.1*kfilt^(-1/3)); % match prefac !!
              iF = Pt*(diag(Fk.^-1)*PW); % note no high freq cutoff,  = A^-1
              F = Pt*(diag(Fk)*PW); Ff = F*f; % F = A, kill high freq by design
              u = iF * (f./xn); Du = D*u; % get neu bdry data from efunc
              A = ixnrm(xn.*u)^2; B = ixnrm(xn.*Du)^2; % get ODE consts
              C = 2*real(ip(Ff,xt.*Du)); E = ixnrm(xn.*Ff)^2;
              khat = evp.solvebetaodeneu(k, be, A, B, C, E); % frozen A,B,C,E
              
            end
            if wantmodes                     % recon boundary efunc...
              if strcmp(o.fhat, 'f'), fhat = f;
              else         % use 1st or 2nd deriv of f
                if ~exist('Df','var'), Df = D*f; mf = m.*f; end 
                Vf = xt.*Df; Vff = ixnip(f,Vf); dfdk = (Vf + mf + Vff*f)/k;
                fhat = f + (khat-k)*dfdk; % 1st, uses best avail khat. Normed
                if strcmp(o.fhat, 's')       % use 2nd-deriv of f:
                % Alex formula, leading order bits only, no c parallel cmpt...
                %  ddfdkk = xn.*xn.*f + (xt.*(D*(Vf+mf)) + m.*(Vf+mf) + ...
                %                        xn.*(D*(xn.*(Df))))/k^2;
                % Andrew's version which drops the negligible xn(xn)'f term...
                  ddfdkk = xn.*xn.*(D*Df/k^2 + f) + ...        % todo: FFTs!
                            (xt.*(D*(Vf+mf)) + m.*(Vf+mf))/k^2;
                  fhat = fhat + (khat-k)^2*ddfdkk/2; % note sign
                end
              end
              fhat = fhat / ixnrm(fhat);   % normalize (avoids needing c)
              nd = sqrt(2) * khat * fhat ./ xn; % normal-deriv function 
            else nd = []; end
            if khat>=kbot & khat<ktop   % if in window, keep it
              kl = [kl; khat]; imbel = [imbel; imag(bes(i))]; % keep Im(beta)
              ndl = [ndl nd]; kol = [kol; k];
            end
          end
          if v, fprintf('found %d eigvals in [%g,%g]\n',numel(kl),kbot,ktop);end
          [dummy i] = sort(kl); % reorder all the l-arrays when append to lists
          kj = [kj; kl(i)]; err.imbej = [err.imbej; imbel(i)];
          err.koj = [err.koj; kol(i)];
          if wantmodes, p.ndj = [p.ndj ndl(:,i)]; end
          kbot = ktop;                        % increment k and repeat
          n = n+1;  % increment data obj counter
        end            % ==== end window loop
        p.kj = kj; p.err = err; p.coj = [];
        
      elseif strcmp(meth,'ms') % ........... SVD iter min search method........
        io = o; io.xtol = o.tol;
        dE = 4*pi/d1.area;  % Weyl mean level spacing in E=k^2
        ppls = 5;          % user param: # grid points per mean levelspacing
        if ~isfield(io, 'maxslope'), io.maxslope = 1.0; end  % not sure why this is good, but O(1)
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
      % OBSOLETE - for Neumann tensions see ./gsvdtension.m
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
    
    function A = fillfredholmop(p, k) % .... helper for solvespectrum
    % FILLFREDHOLMOP - fill l^2 normed (sqrtwei) layer-potential matrix at k
    %
    % A = FILLFREDHOLMOP(p, k) fills the matrix A = (I - 2D) for Dirichlet BCs
    %    or -(I + 2D^T) for Neumann BCs, at wavenumber k, for evp object p.
    %    The appropriate layerpotential must have been set up in the evp object.
    %
    % This works for Dirichlet (bc='D', lp='d') or Neumann (bc='N', lp='s').
    % It is a helper routine for solvespectrum meth='fd' and 'ms', and
    % for solvemodecoeffs meth='ms'.
    %
    % Issues: not fully documented yet.
      if numel(k)~=1, error('k must be a single number!'); end
      p.setoverallwavenumber(k);
      A = p.fillbcmatrix;   % fill matrix mapping density values to field values
      N = size(A,1); if size(A,2)~=N, error('A must be square!'); end
      % rescale so l^2 norms represent L^2 norms, as Jim Bremer also does...
      A = A .* repmat(1./p.sqrtwei, [N 1]); % note sqrtwei already left-applied 
      A = -2*A;     % for Dirichlet, -2 turns (D-1/2) into (1-2D) for det
    end
    
    function [coj ndj minsigj] = solvemodescoeffs(p, meth, o) % all evals in pr
    % SOLVEMODESCOEFFS - find modes and n-derivs for all eigenwavenumbers in p
    %
    % [coj ndj minsigj] = SOLVEMODESCOEFFS(p, meth) takes eigenwavenumbers kj in
    %   given p evp object and, starting with only this information, computes an
    %   eigenmode representation (coefficients and boundary data) at each.
    %   Degeneracies are identified (see o.degtol below), and for them, a set of
    %   linearly-independent modes spanning the eigenspace is returned.
    %   Solution method is as in SOLVEMODECOEFFS. Cost is O(N^3) per mode.
    %   An appropriate layer-potential basis must already exist in p.
    %
    %   Outputs:
    %   coj - stack of column vecs of coefficients of each mode,
    %   ndj - stack of col vecs giving boundary function data for each mode,
    %   minsigj - minimum singular values for each kj: small values verify kj
    %             as a good approximate eigenfrequency.
    %   Note: this data is also saved in p object, overwriting p.coj, p.ndj and
    %         p.err.minsigj.
    %
    % [...] = SOLVEMODESCOEFFS(p, meth, opts) controls options such as:
    %   opts.verb = 0,1,... : verbosity (0=silent, 1=diagnostic text, ...)
    %   opts.degtol : absolute tolerance on closeness of kj to count as
    %                 degenerate (default 1e-12). Choose to match desired
    %                 tolerance on kj.
    %
    % See also: EVP.SOLVESPECTRUM, EVP.SOLVEMODECOEFFS
      if nargin<2, meth = []; end
      if nargin<3, o = []; end
      if ~isfield(o, 'verb'), o.verb = 1; end; v = o.verb;  % default verbosity
      if ~isfield(o, 'degtol'), o.degtol = 1e-12; end;  % default tol
      ne = numel(p.kj);
      if ne<1, warning('no eigenwavenumbers kj in evp object!'); end
      j=1; p.coj = []; p.ndj = []; p.err.minsigj = [];  % erase existing data
      while j<=ne
        dim = numel(find(abs(p.kj(j:end)-p.kj(j))<o.degtol)); % degeneracy
        js = j:j+dim-1;  % list of indices in this subspace
        [co nd e] = solvemodecoeffs(p, p.kj(js), meth, o);
        p.coj = [p.coj co]; p.ndj = [p.ndj nd];  % stack (todo: realloc?)
        p.err.minsigj = [p.err.minsigj; e];
        if v, k = p.kj(js(1));  % text output (k thereof)
          if dim==1, fprintf('   mode #%d      \tk=%.16g: min sing val=%.3g\n',j,k,e);
          else, fprintf('   modes #%d-%d   \tk=%.16g: sing vals in [%.3g,%.3g]\n',min(js),max(js),k,min(e),max(e));
          end
        end
        j = j+dim;
      end
      if nargout>=1, coj = p.coj; end        % outputs if needed...
      if nargout>=2, ndj = p.ndj; end
      if nargout>=3, minsigj = p.err.minsigj; end
    end
    
    function [co nd e] = solvemodecoeffs(p, k, meth, o) % ........ single efunc
    % SOLVEMODECOEFFS - compute mode or eigenspace for single wavenumber; O(N^3)
    %
    % Compute an eigenmode representation (coefficients and boundary data) at
    %   a single eigenfrequency k, or eigenspace at at degenerate cluster of k,
    %   using the minimum singular value(s) of 1/2-D (for Dirichlet BC case)
    %   or 1/2+D^* (for Neumann case).
    %   Desired eigenspace dimension (dim) is indicated by passing in a list of
    %   (repeated) values for k.
    %   An appropriate layer-potential basis must already exist in p.
    %
    % Outputs: co - (N-by-dim) coeff data, ie columns of density functions
    %                          (single-layer for Neu, double-layer for Dir).
    %          nd - (N-by-dim) columns of mode boundary data (normal-derivative
    %                          for Dir; boundary values for Neu).
    %           e - (dim-by-1) singular value(s) of Fredholm operator - should
    %               be small if k is a good approximate eigenwavenumber.
    %
    % [...] = solvemodecoeffs(p, k, meth, opts) sets various options, including:
    %    opts.meth = 'ms' (default; the only implemented method) use min sing
    %              values of Fredholm operator (I - 2D) for Dirichlet BCs, or
    %              (I + 2D^T) for Neumann BCs. Operator chosen by fillfredholmop
    %    opts.iter = 0 use matlab's full svd to get the lowest singular values
    %                  and vectors,
    %                1 use custom iterative method (default), 10x faster.
    %                  (For degenerate cases, reverts to iter=0 method.)
    % Notes/issues:
    % * only Dirichlet BC is normalized correctly, since this uses Rellich
    %   formula. This could be done for general BCs, would need specrowdiff to
    %   compute tangential derivs. However, we don't know how to normalize
    %   densities anyway, and SLP is used for Neumann eval since better.
    % * Braxton Osting needed domains w/ multiple segments, so that is done.
    %
    % See also: EVP.FILLFREDHOLMOP

    if nargin<3 | isempty(meth), meth='ms'; end   % default method
      if nargin<4, o = []; end
      if ~isfield(o, 'iter'), o.iter = 1; end   % default iterative
      
      if strcmp(meth,'ms') %....... min sing val(s) method
        
        dim = numel(k); if dim>1, k=mean(k); end
        A = p.fillfredholmop(k);      % eg (I - 2D) in Dirichlet, (I + 2D^T) Neu
        if o.iter && dim==1
          [u e v info] = utils.minsingvalvecs(A);  % iter O(N^3) small const
          if info.flag
            warning('info.flag failed from utils.minsingvalvec, probably due to near-degeneracy: switching to opts.iter=0 instead for this eigenspace');
            [U S V] = svd(A); u = U(:,end-dim+1:end); v = V(:,end-dim+1:end);
            s = diag(S); e = s(end-dim+1:end);   % O(N^3) 10x slower method
          end
        else                                       % O(N^3) 10x slower method
          [U S V] = svd(A); u = U(:,end-dim+1:end); v = V(:,end-dim+1:end);
          s = diag(S); e = s(end-dim+1:end);
        end
        % co = density, converted from l^2 to value vector (unnormalized)
        co = v .* repmat(1./p.sqrtwei.', [1 dim]);
        % L sv's are R sv's w/ D^*, so u gives boundary function...
        x = vertcat(p.segs.x); nx = vertcat(p.segs.nx); % allow multiple segs
        xdn = real(conj(x).*nx); % Rellich bdry wei aka Morawetz multiplier
        for i=1:dim
          u(:,i) = u(:,i)/max(u(:,1)); % make Re-valued (Im part is error est)
          if p.segs(1).a~=0   % assume Dirichlet BC (everywhere)
            rellichint = sum(xdn.*abs(u(:,i)).^2);      % Rellich integral
          else                % Neumann BC (assumed everywhere)
            rellichint = 1;  % TODO: build D matrix on multi segs...
            % use Barnett generalized inner product identity...
            %Du = p.sqrtwei.' .* (D*(u(:,i)./p.sqrtwei.')); %tang deriv * sqrtwei
            %rellichint = sum(xdn.*(k(1)^2*abs(u(:,i)).^2 - abs(Du).^2));
          end
          u(:,i) = (sqrt(2) * k(1)) * u(:,i) / sqrt(rellichint);
        end
        nd = u .* repmat(1./p.sqrtwei.', [1 dim]); % convert l^2 to values
      
      else
        error('unknown method');
      end
    end
    
    function [e] = modeserrors(p, q, o) % ... L2 errors of modes vs ref set
    % MODESERRORS - L2 bdry errors of set of modes relative to a reference set
    %
    % e = MODESERRORS(p, q) where p and q are evp objects computes
    %   angles between boundary functions of modes in p against corresponding
    %   modes in q. p and q must contain the same geometry and set of modes,
    %   otherwise the results are meaningless.
    %   In case of degeneracies, the principal subspace angle between
    %   eigenspaces is used, and this one value is written out to all entries
    %   in the cluster. See options for choice of L2 inner product weight.
    %
    % e = MODESERRORS(p, q, opts) controls options such as:
    %   opts.degtol : absolute tolerance on closeness of p.kj to count as
    %                 degenerate (default 1e-12). Choose to match desired
    %                 eigenwavenumber tolerance.
    %  opts.wei : 0 uses plain L2(bdry) applied to d_n phi (default)
    %             1 uses Rellich weight (x.n), equivalent to <f,f> inner prod
    %               from Barnett-Hassell scaling paper (must be star-shaped).
      if nargin<3, o = []; end
      if ~isfield(o, 'degtol'), o.degtol = 1e-12; end;
      if ~isfield(o, 'wei'), o.wei = 0; end;
      ne = numel(p.kj);
      if numel(q.kj)~=ne, error('different numbers of modes, stuck!'); end
      if size(q.ndj)~=size(p.ndj), error('ndj different sizes, stuck!'); end
      w = p.sqrtwei.';   % vanilla L2 weighting for derivative values in ndj
      if o.wei==1
        x = vertcat(p.segs.x); nx = vertcat(p.segs.nx); % allow multiple segs
        xdn = real(conj(x).*nx); w = w.*sqrt(xdn); % l2 wei for phi_n IP
      end
      j=1; e = nan(size(p.kj));
      while j<=ne
        dim = numel(find(abs(p.kj(j:end)-p.kj(j))<o.degtol)); % degeneracy
        js = j:j+dim-1;  % list of indices in this subspace
        if max(q.kj(js))-min(q.kj(js))>o.degtol
          warning(sprintf('eigenwavenumbers #%d-%d not degenerate (degtol=%g) in evp object q', min(js), max(js),o.degtol));
        end
        e(js) = subspace(p.ndj(:,js).*repmat(w,[1 dim]), q.ndj(:,js).*repmat(w,[1 dim]));
        j = j+dim;
      end
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
    %   opts.eval - 'basis' uses existing basis p.coj to eval (default for Neu)
    %               'grf' uses Green's rep. formula on p.ndj (default for Dir)
    %   opts.kwin - only computes modes in given wavenumber window
    %   opts.inds - only computes modes in given (unordered) index set
    %   opts.col - 'jet' (default, matlab colorscale), 'bw' (phi^2 in grey)
    %   opts.nac - override number of subplots across
    %
    % Notes/issues:
    %  * Normalization hasn't been considered much. For GRF case, they are
    %    correctly L2-normalized over the domain. Phase removal is a bad hack!
    %  * No account taken of non-simply connected domains (needs GRF for
    %    this to be set up...)
    %  * mixed D-N domains, Robin ... ?
    %  * Close-evaluation (J-expansion projection) for layer-potentials...?
    %  * There is a hack to use an interior pt for phase normalization - fix.
    %
    % See also: EVP, PROBLEM.SHOWSOLUTION
      if nargin<2, o = []; end
      neubc = (p.segs(1).a==0);  % assume Neumann BCs everywhere
      if ~isfield(o, 'eval')
        if neubc, o.eval='basis'; else, o.eval='grf'; end
      end
      grf = strcmp(o.eval,'grf');
      if ~isfield(o, 'kwin'), o.kwin=p.kwin; end
      if ~isfield(o, 'inds'), o.inds=1:numel(p.kj); end
      if ~isfield(o, 'bdry'), o.bdry = 0; end  % unused: future expansion
      if ~isfield(o, 'col'), o.col = 'jet'; end
      wantdata = nargout>0;
      inlist = 0*p.kj; inlist(o.inds) = 1;         % true if in index list
      js = find(p.kj>o.kwin(1) & p.kj<o.kwin(2) & inlist); % indices to plot
      if isempty(js)
        warning('no eigenvalues found in evp, or list or window!'); end      
      
      if grf
        if isempty(p.ndj), error('GRF method: no ndj bdry func info found!');end
        d = utils.copy(p.doms);  % set up for GRF evaluations in loop...
        d.clearbases; if neubc, d.addlayerpot(d.seg, 'd'); else
          d.addlayerpot(d.seg, 's'); end  % correct GRF type for BCs
        pe = bvp(d); pe.setupbasisdofs; % temporary new problem instance
      else
        if isempty(p.coj), error('basis method: no coj coeffs info found!');end
        pe = p;                  % just use existing problem
      end
      n = numel(js);
      figure; nac = ceil(sqrt(n)); if isfield(o,'nac'), nac = o.nac; end
      ndn = ceil(n/nac); % # subplots across & down

      for i=1:n, j = js(i);      % ---- loop over selected eigenvalues
        pe.setoverallwavenumber(p.kj(j));                         % get this k
        if grf, pe.co = p.ndj(:,j); else, pe.co = p.coj(:,j); end % get coeffs
        [u gx gy di] = pe.gridsolution(o); % either conventional or GRF evalu
        %if ~grf  % rotate phase of u (est at int pt) so that essentially real
        %  uint = u(floor(numel(gx)/2),floor(numel(gy)/2)) % HACK guess int pt!
        %  u = u * conj(uint/abs(uint));  % unphase u
        %end
        u = real(u);                     % saves half the memory (real-valued)
        dxdy = (gx(2)-gx(1))*(gy(2)-gy(1));
        if ~grf | neubc       % we don't know how to normalize the basis eval
          u = u / sqrt(sum(dxdy*abs(u(~isnan(u))).^2)); % grid normalize crudely
        end
        if wantdata
          if i==1, uj = nan(size(u,1), size(u,2), n); end % allocate output
          uj(:,:,i) = u;                                  % copy data to array
        end
        tsubplot(ndn, nac, i);
        if strcmp(o.col,'jet')
          imagesc(gx, gy, u, 'alphadata', ~isnan(di)); colormap(jet(256));
 % contourf(gx, gy, u, [-1:0.2:1]*abs(u(ceil(end/2),ceil(end/2)))); colormap(jet(256)); hold on; p.segs.plot; % debug with known vertical scale
          caxis(3.5*[-1 1]/sqrt(p.doms.area)); % rescale values based on area
        elseif strcmp(o.col,'bw')
          imagesc(gx, gy, u.^2, 'alphadata', ~isnan(di)); colormap(1-gray(256));
          caxis([0 4]/p.doms.area); % rescale values based on area
        end
        axis equal tight; axis off;
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
    [d V] = filteredDtNspectrum(p, kstar, opts)
   end % methods
  
  % --------------------------------------------------------------------
  methods(Static)    % these don't need EVP obj to exist to call them...
    [k_weyl] = weylcountcheck(k_lo, k, perim, area, dkfig)   % old Weyl routine
    kh = solvebetaode(kstar, beo, A, B, C, Ap, Bp, Cp)  % NtD scaling helper
    kh = solvebetaodeneu(kstar, beo, A, B, C, E)  % DtN scaling helper
    [xm ym info] = gridminfit(f, g, o)   % gridded vectorial minimizer
    [xm fm info] = iterparabolafit(f, x, y, opts)   % minimization helper
    [A,B,C] = para_fit(e, f)        % parabolic helper
    [t V F G A B] = tensionsq(d, E, opts)  % domain eigenvalue solver helper
    Fh = spectralfiltermatrix(s,k,opts)    % Neumann inclusion helper for ninc
    [t V] = gsvdtension(d, k, opts)    % Neumann eigval finding helper for ninc
    [G An] = intnormmatrix(s, k, A, Dx, Dy)   % helper for gsvdtension
  end
end
