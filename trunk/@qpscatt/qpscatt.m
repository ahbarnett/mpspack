% QPSCATT - define quasi-periodic (grating) frequency-domain scattering problem
%
%  pr = qpscatt(airdoms, doms, d) creates a scattering problem object pr
%   as with scattering, except a quasi-periodic version with x-periodicity of d.
%   The solution method will be via FTyLPs, with possibly shifted Sommerfeld
%   contour and pole correction.
%   airdoms must be only the single exterior domain of the obstacle scattering
%   problem.
%
%   Note: currently only layer-potential bases are allowed in the exterior
%   domain of the obstacle (adding the evalfty method to other bases would be
%   needed).
%
%  pr = qpscatt(airdoms, doms, d, opts) sets certain algorithm parameters:
%    opts.M = # degrees of freedom in each FTyLP
%    opts.nei = 0,1,2... how many direct sums to include on each side
%    opts.buf = 0,1,2... how wide the qpstrip domain is including buffer size
%
% See also PROBLEM, BVP, SCATTERING

% Copyright (C) 2010 Alex Barnett

classdef qpscatt < scattering & handle
  properties
    a                                       % alpha for problem periodicity
    t                                       % the qpstrip object (has own alpha)
    d                                       % problem periodicity (not qpstrip)
    extdom                                  % the single exterior domain to obst
    nei                                     % 0,1,...: 1,3,... direct images
    buf                                     % 0,1,...: qpstrip 1,3,... periods
    M                                       % # nodes FTyLP Sommerfeld contour
    kpx, kpy                                % k-plane poles (k_x,k_y) (QP waves)
    kpn                                     % Bragg orders of these poles
    ikpfix                                  % indices of k_y for poles to cross
    B, T                                    % segments T,B for Bragg ampl fix
  end
  
  methods % -------------------------------- methods particular to qpscatt
    function pr = qpscatt(airdoms, doms, d, o) % ................... constructor
      if nargin==0, airdoms = []; doms = []; end   % empty constructor
      pr = pr@scattering(airdoms, doms);           % cannot be conditionalized
      if ~isempty([airdoms doms])           
        if numel(airdoms)~=1, error('airdoms must be a single domain');end
        if ~airdoms.exterior, error('airdoms must be an exterior domain');end
        pr.d = d; pr.extdom = airdoms; pr.a = 1.0;   % default
        if nargin<4, o = []; end            % process method options
        if isfield(o,'M'), pr.M = o.M; else pr.M = 120; end  % defaults
        if isfield(o,'nei'), pr.nei = o.nei; else pr.nei = 2; end
        if isfield(o,'M'), pr.buf = o.buf; else pr.buf = 0; end
        if pr.nei<pr.buf, warning('nei<buf: probably will not work!'); end
        % create a qpstrip... strip vector t.e may be larger than periodicity d
        pr.t = qpstrip(d*(1+2*pr.buf), 1.0);            % dummy/default omega 
        % choose a periodizing basis for the qpstrip...
        pr.t.addqpftylayerpots(struct('M',pr.M,'nearsing',3)); % NB default!
        % create sensible B, T segments...
        Mt = 20; safedist = 0.5;          % Mt, safedist could dep on omega?
        yB = min(imag(airdoms.x)) - safedist;
        yT = max(imag(airdoms.x)) + safedist;
        pr.B = segment(Mt, 1i*yB+[-d/2,d/2], 'p');
        pr.T = segment(Mt, 1i*yT+[d/2,-d/2], 'p');
        pr.T.x = pr.T.x + pr.T.w(1)/2;
        pr.B.x = pr.B.x - pr.B.w(1)/2; % safely away from walls 1/2 a gridpoint
      end
    end

    function setoverallwavenumber(pr, k) % ................. overloads @problem
    % SETOVERALLWAVENUMBER - set problem k in each domain using refractive index
    %
    %  setoverallwavenumber(pr, omega) propagates overall wavenumber omega to
    %   each domain. If a domain has refractive index n, then
    %   its wavenumber will become n.omega. The qpstrip also is set to omega.
      setoverallwavenumber@problem(pr, k);
      pr.t.k = k;                            % set the qpstrip
    end
    
    function setincidentwave(pr, t)
    % SETINCIDENTWAVE - set inc. plane wave direction, Bloch, k-fix, QP bases
    %
    %  setincidentwave(pr, t) sets up the Bloch alpha and BCs required in
    %   qpscatt problem pr, for planewave at angle t in [-pi,0].
    %   Wavenumber must already be set in the problem. FTyLPs are adjusted
    %   based on nearsing, and Wood's correction bases are chosen and set up.
    %
    %  Careful: calling this routine overwrites all inhomogeneity functions or
    %   data f, g stored on any of the problem's segments. However it preserves
    %   existing a, b BC or matching coeffs (which must be set up on entry).
      if nargin==2
        if isempty(pr.k), error('qpscatt problem needs wavenumber set'); end
        t = mod(t,2*pi);
        if t<pi, error('upwards incident angle (use downwards)!');end
        setincidentwave@scattering(pr, t); % call superclass method
        om = pr.k; kvec = om*exp(1i*t); pr.a = exp(1i*real(conj(kvec) * pr.d));
        pr.t.setbloch(exp(1i*real(conj(kvec) * pr.t.e))); % set Bloch alpha
        % set the k-poles list in pr object...
        n = 2 * ceil(om*pr.d/2/pi) + 5; kpn = -n:n; % hack # diffraction orders?
        kpx = om*cos(t)+kpn*2*pi/pr.d; kpy = sqrt(om^2-kpx.^2); % kp y-compt
        [dum i] = sort(abs(kpy),'ascend');       % note +ve sqrt key here
        pr.kpx = kpx(i); pr.kpy = kpy(i); pr.kpn = kpn(i); % sorted lists
        pr.ikpfix = []; shift = 0; toonear = 1.4; % dist scale for Wood anom.
        k1 = abs(pr.kpy(1)); k2 = abs(pr.kpy(2)); % two poles closest to origin
        if k1<toonear                % polefix decision regions for (k1,k2)
          pr.ikpfix = 1;
          if k2<3*toonear
            shift = max(2*toonear, k2+toonear); pr.ikpfix(2) = 2;
          else, shift = 2*toonear; end
        end
        kde = min(abs(shift-abs(pr.kpy)))/sqrt(2); % estimate closest approach
        pr.t.bas = pr.t.bas(1:2); % reset QP FTyLP bases, then add new ones...
        for i=1:2, pr.t.bas{i}.requadrature(pr.t.bas{i}.lp.N, ...
                      struct('omega',om, 'shift', shift, 'nearsing', 2.0*kde));
        end
        N = pr.t.bas{1}.lp.Nf; % now measure dist(shifted contour, poles)...
        kd = min(min(abs(repmat(pr.kpy.', [1 N]) - ...
                         repmat(pr.t.bas{1}.lp.kj,[numel(pr.kpy) 1]))));
        fprintf('# poles fixed = %d, shift = %.3g, est dist = %.3g, dist = %.3g\n', numel(pr.ikpfix), shift, kde, kd)
        if kd<toonear, warning('dist(shifted contour, poles) = %.3g\n', kd);end
        if numel(pr.ikpfix)>0 % add any Wood's anomaly fix bases...
          pr.t.addepwbasis(1, struct('real',0,'half',0));
          pr.t.bas{3}.dirs=(pr.kpx(pr.ikpfix(1))+1i*pr.kpy(pr.ikpfix(1)))/om;
        end
        if numel(pr.ikpfix)>1
          pr.t.addepwbasis(1, struct('real',0,'half',0));
          pr.t.bas{4}.dirs=(pr.kpx(pr.ikpfix(2))+1i*pr.kpy(pr.ikpfix(2)))/om;
        end
      end
    end
    
    function showkyplane(pr)
    % SHOWKYPLANE - plot Sommerfeld quadrature nodes, k_y poles, +-om, etc
      figure; plot(pr.t.bas{1}.lp.kj, '+'); axis tight equal; hold on;
      plot([pr.kpy -pr.kpy], 'rx');
      i = pr.ikpfix; plot(real(pr.kpy(i)),imag(pr.kpy(i)),'go');
      plot([-pr.k pr.k], [0 0], 'k*');
      title('complex k (ie k_y) plane: contour, poles'); 
    end
    
    function showbragg(pr)
    % SHOWBRAGG - plot upwards Bragg directions and theta_inc, etc
      plot([0 cos(pr.incang)],[0 sin(pr.incang)], 'k-', 'linewidth', 3);
      i = find(imag(pr.kpy)==0);    % indices of propagating Bragg modes
      o = 0*pr.kpx(i);              % row of zeros for starting points
      hold on; plot([o; pr.kpx(i)/pr.k],[o; pr.kpy(i)/pr.k],'m-','linewidth',3);
    end
    
   
    function rhs = fillrighthandside(pr)
    % FILLRIGHTHANDSIDE - overloads BVP routine, setting mismatch and discrep
    %
    % Note: incident field is QP so leads to zero discrepancy.
      rhs = [];                 % get ready to stack stuff onto it as a big col
      obstrhs = fillrighthandside@bvp(pr); % obstacle mismatch part of RHS
      Nq = pr.t.setupbasisdofs;
      rhs = [obstrhs; zeros(Nq, 1)];       % QP block square, so use # QP dofs
      if nargout==0, pr.rhs = rhs; end     % this only stores pr.rhs if no outp
    end

    function A = fillbcmatrix(pr, opts)
    % FILLBCMATRIX - fills [A B;C Q] for QP scatt prob: dofs->mismatch/discrep
    %
    % Issues/Notes:
    %  * obst bases affecting extdom (airdoms) can only be LPs (since need to
    %    use evalfty). Would need evalfty methods for other exterior basis sets.
    %  * everything dense. Need replace Q by diagonal multiplication
    %    method, and make A a method not a matrix for FMM case.
    %    Better, make A a struct with mult method, B, C, and diag(Q) in it.
      bob = pr;                        % following fillbcmatrix@problem
      N = bob.setupbasisdofs;          % sets up obst dofs in pr
      M = numel(pr.sqrtwei);
      t = pr.t;                        % the qpstrip domain
      Nq = pr.t.setupbasisdofs;        % qpstrip domain's QP bases
      Mq = pr.t.bas{1}.Nf;             % # FTyLP quadr pts, 1/2 the discrep dofs
      nb = numel(pr.t.bas);            % how many QP bases in t
      nei = pr.nei; buf = pr.buf; d = pr.d; a = pr.a; % get Bloch alpha
      
      % fill Q...
      Q = t.evalbasesdiscrep();        % is block-diag - need not really fill!
      
      % fill A: self, then directly sum neighbor images contribs...
      A = fillbcmatrix@problem(pr);    % obstdofs->mismatch block
      for n=[-nei:-1 1:nei]
        Ac = fillbcmatrix@problem(pr, struct('trans',-d*n, 'doms',pr.extdom));
        A = A + a^n * Ac;
      end
      % fill B: effect of FTyLPs on the airdom-touching segs in pr...
      B = zeros(M, Nq); m = 0;                                  % B row counter
      for s=pr.segs           % loop over segs -------------------
        if s.bcside==0                       % matching, may be dielectric
          if s.dom{1}.isair~=s.dom{2}.isair    % air-nonair (ext-nonext) junct
            ms = m + (1:2*size(s.x,1));     % 2M colloc indices for block row
            % extract relevant coeffs for whichever side have inc field on it:
            if s.dom{1}.isair, sa = s.a(1); sb = s.b(1); % affects + side only
            else sa = s.a(2); sb = s.b(2); end % affects - side only
            for i=1:nb, b = t.bas{i}; ns = t.basnoff(i)+(1:b.Nf);
              [Bb Bbn] = b.eval(s);          % val and nderiv both always need
              if isa(b, 'ftylayerpot') % explicitly make copies on R
                [Bbt Bbtn] = b.eval(s.translate(-t.e));
                Bb = Bb + a*Bbt; Bbn = Bbn + a*Bbtn;
              end
              % TODO: could speed up by testing for nonzero sa, sb...
              B(ms,ns) = repmat(pr.sqrtwei(ms).', [1 b.Nf]) .* ...
                  [sa * Bb; sb * Bbn];  % stack value then deriv contribs
            end
          end
        elseif s.bcside==1 | s.bcside==-1    % BC
          ind = (1-s.bcside)/2+1;   % index 1 or 2 for which side the BC on
          if s.dom{ind}.isair       % relevant domain, air-to-metallic boundary
            ms = m+(1:size(s.x,1)); % colloc indices for this block row
            for i=1:nb, b = t.bas{i}; ns = t.basnoff(i)+(1:b.Nf);
              if s.b==0             % Dirichlet only
                Bblock = s.a * b.eval(s);
                if isa(b, 'ftylayerpot') % explicitly make copies on R
                  Bblock = Bblock + s.a*a*b.eval(s.translate(-t.e));
                end
              else                  % Robin, Neumann
                [Bb Bbn] = b.eval(s);
                if isa(b, 'ftylayerpot') % explicitly make copies on R
                  [Bbt Bbtn] = b.eval(s.translate(-t.e));
                  Bb = Bb + a*Bbt; Bbn = Bbn + a*Bbtn;
                end
                Bblock = s.a * Bb + s.b * Bbn;  % value and deriv contribs
              end
              B(ms,ns) = repmat(pr.sqrtwei(ms).', [1 b.Nf]) .* Bblock;
            end
          end
        end
        m = ms(end);                    % update counter in either case
      end % ------------------ segs loop
      
      % fill C: discrep of airdoms-affecting obst bases (must be LPs!) in pr...
      C = zeros(2*Mq, N); f = utils.copy(t.bas{1}.lp); % dummy bas, transl evals
      o = []; o.dom = pr.extdom; % must be a single exterior domain
      for i=1:numel(bob.bas)           % all bases in problem or basis obj
        b = bob.bas{i}; ns = bob.basnoff(i)+(1:b.Nf); % dof indices for this bas
        if utils.isin(pr.extdom, b.doms) % restrict to bases affecting airdoms
          for n=nei-2*buf:nei   % general buf case (nei>=buf)
            f.orig = t.Lo-n*d; [Cj Cxj] = b.evalfty(f, o);
            C(:,ns) = C(:,ns) + a^n*[Cj;Cxj];
            f.orig = t.Lo-(-1-2*buf-n)*d; [Cj Cxj] = b.evalfty(f, o);
            C(:,ns) = C(:,ns) - a^(-1-2*buf-n)*[Cj;Cxj];
          end  % careful, since L's phase = 1 always
        end
      end  % bob basis loop
      
      pr.A = [A B; C Q];       % stack the full E matrix (it's called A)
      
      if ~isempty(pr.ikpfix)   % augment matrix w/ extra rows for Wood's fix ...
        [RO RI] = pr.fillbraggamplrow(pr.ikpfix);   % get Bragg ampl matrix
        RI = RI.*repmat(1./sqrt(sum(RI.^2, 2)), [1 size(RI,2)]); % row-normalize
        pr.A = [pr.A; RI];     % append block row to E
      end  
    end % fillbcmatrix
    
    function [RO RI] = fillbraggamplrow(pr, i, o)
    % FILLBRAGGAMPLROW - row(s) of full E matrix measuring Bragg amplitude(s)
    %
    % [RO RI] = FILLBRAGGAMPLROW(pr, i) where pr is a qpscatt problem, returns
    %  the row(s) of the full (obst + QP bases) E matrix which maps to the
    %  outgoing and incoming Bragg amplitude(s) indexed by i, on segment pr.S.
    %  The indexing of i is the same as that of pr.kpx, pr.kpy, and pr.kpn.
    %
    % [RO RI] = FILLBRAGGAMPLROW(pr, i, opts) allows options such as,
    %   opts.side = 'T', 'B' (default): choose if meas on top or bottom seg
    %   opts.test: if present, self-test using incident wave
    %
    % Notes: adapted from code polesftylp.m. k-scaling removed, to meas coeffs
      if nargin<3, o = []; end
      if ~isfield(o, 'side'), o.side = 'B'; end
      if o.side=='B', S = pr.B; elseif o.side=='T', S = pr.T;
      else error('invalid opts.side'); end
      [AS ASn] = pr.evalbases(S); % evals all bases (obst + QP) on S pointset
      if isfield(o,'test'), AS = pr.ui(S.x); ASn = pr.uix(S.x).*real(S.nx) + pr.uiy(S.x).*imag(S.nx); end % test w/ inc wave, should give abs(ampl)=1
      % project onto x-Fourier modes (relative to origin)...
      ni = numel(i);
      F = repmat(S.w, [ni 1]).*exp(-1i*pr.kpx(i).'*real(S.x.')); % outerprod
      FS = F*AS; FSn = F*ASn;                       % each numel(i)-by-(#dofs)
      % use little 2x2 on each Bragg index i to map to up & down-going ampls...
      %pref = 0.3;                                  % by-hand row scale factor
      RO = zeros(numel(i), size(AS,2)); RI = RO;    % allocate matrices
      for j=1:numel(i), k = pr.kpy(i(j));           % loop over kpy indices...
        if k==0, ab2ce = eye(2)/pr.d;               % the 2x2: better k->0 idea?
        else ab2ce = [1 1/(1i*k);1 -1/(1i*k)]/2/pr.d; end
        Rce = ab2ce * [FS(j,:); FSn(j,:)];          % NB row j not row i
        RO(j,:) = Rce(1,:); RI(j,:) = Rce(2,:);     % copy into output arrays
      end
    end
    
    function [bo to bi ti] = braggampl(pr, i, o)
    % BRAGGAMPL - compute (subset of) scatt field Bragg amplitudes via B,T segs
    %
    % [bo to bi ti] = BRAGGAMPL(pr) returns Bragg amplitudes for all orders
    %   present in pr.kpn Bragg order list. They are computed using matrices
    %   evaluated on (asssumed sensible) existing pr.B and pr.T segments.
    %   pr is the qpscatt problem object, which must have coefficients pr.co.
    %
    % Note: this is for scattered field not the full field.
    %
    % See also: FILLBRAGGAMPLROW
      if nargin<3, o = []; end
      if isempty(pr.co), error('coefficient vector is empty!'); end
      if nargin<2 || isempty(i), i = 1:numel(pr.kpn); end % default is all i
      [RO RI] = pr.fillbraggamplrow(i, o);   % get B (bottom) Bragg ampl matrix
      co = pr.co; if isfield(o, 'test'), co = 1; end % dummy for test
      bo = RO*co; bi = RI*co;
      o.side = 'T'; [RO RI] = pr.fillbraggamplrow(i, o);   % T (top)
      to = RO*co; ti = RI*co;
    end
    
    function [u d n] = braggpowerfracs(pr, o)
    % BRAGGPOWERFRACS - compute propagating Bragg scattered power fractions
    %
    % [u d n] = BRAGGPOWERFRACSx(pr) returns upwards and downwards scattered
    %   intensities u, d (for unit incident beam intensity) for periodic
    %   scattering problem pr, given coefficients pr.co. n is respective Bragg
    %   orders list. (u,d,n are column vectors)
    %
    %   opts.test: if present, self-test using incident wave; gives u=d=0
    %   opts.table: if present, print a table showing coeffs
    %
    % See also: BRAGGAMPL
      if nargin<2, o = []; end
      i = find(imag(pr.kpy)==0);        % indices of propagating orders
      if isfield(o,'test'), i = find(pr.kpn==0); end % fake for single inc wave
      n = pr.kpn(i)';                   % col vec
      % flux angle factors propto k_y (rescaled to 1 if same as inc)...
      angfacs = pr.kpy(i); angfacs = angfacs/angfacs(find(pr.kpn(i)==0));
      [bo to bi ti] = pr.braggampl(i, o);
      j0 = find(pr.kpn(i)==0);          % index within i of incident (0) order
      bo(j0) = bo(j0) + pr.ui(1i*imag(pr.B.x(1))); % add inc PW ampl
      ti(j0) = ti(j0) + pr.ui(1i*imag(pr.T.x(1))); % " (irrelevant, ti not used)
      u = to.*conj(to).*angfacs.';
      d = bo.*conj(bo).*angfacs.';
      if isfield(o,'table')           % show results table and flux error
        fprintf('\t\tBragg order\tflux fracs: up\t\t\tdown\n')
        disp([n u d]); fprintf('Bragg tot flux err = %.3g\n',sum([u;d])-1)
      end
    end
    
    function [u di] = pointsolution(pr, p) % ............... overloads @problem
    % POINTSOLUTION - evaluate soln to periodic problem a pointset, given coeffs
    %
    % Overloads problem.pointsolution, but creating dummy problem and dummy
    %  extdom which marries the bases from the obst extdom and the qpstrip,
    %  and also sums in the neighboring direct images of the extdom bases.
    %
    % See also PROBLEM.POINTSOLUTION, QPSCATT
    
      % make dummy problem with new extdom object (whose d.bas are as before)
      pd = utils.copy(pr); pd.A = []; % (a hack!). leave p.co since need it
      for i=1:numel(pd.doms), d = pr.doms(i); if d.isair, 
          d = utils.copy(d); pd.doms(i) = d; pd.extdom = d; end, end
      % now add the QP bases onto those in obst exterior domain...
      pd.extdom.bas = {pd.extdom.bas{:} pd.t.bas{:}}; % marry the bases lists
      pd.bas = {pr.bas{:} pd.t.bas{:}}; % append problem bases ...IN ORDER!
      pd.basnoff = [pr.basnoff pr.N+pr.t.basnoff]; % effective setupbasisdofs
      pd.N = pr.N + pr.t.N;            % finish up appending to pd basis setup
      [u di] = pointsolution@problem(pd, p);     % pass to usual evaluator
      extdomi = find([pr.doms.isair]);           % index of the extdom
      if numel(extdomi)~=1, error('there appears to be >1 exterior domain!');end
      for n=[-pr.nei:-1 1:pr.nei]                % direct sum nei contribs
        uc = pointsolution@problem(pr, pointset(p.x-pr.d*n, p.nx));
        u = u + pr.a^n * uc .* (di==extdomi);    % only affects obst exterior
      end
    end

    function [A Ax Ay] = evalbases(pr, p, opts) % ........ eval obst + QP bases
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
    %   object, jump relations will be taken into account if opts.dom is given.
    %
    % See also: PROBLEM.EVALBASES, QPSCATT.POINTSOLUTION, PROBLEM.DOMAININDICES
      if nargin<3, opts = []; end
      % make dummy problem with new extdom object (whose d.bas are as before)
      pd = utils.copy(pr); pd.A = []; % (a hack!). leave p.co since need it
      for i=1:numel(pd.doms), d = pr.doms(i); if d.isair, 
          d = utils.copy(d); pd.doms(i) = d; pd.extdom = d; end, end
      % now add the QP bases onto those in obst exterior domain...
      pd.extdom.bas = {pd.extdom.bas{:} pd.t.bas{:}}; % marry the bases lists
      pd.bas = {pr.bas{:} pd.t.bas{:}}; % append problem bases ...IN ORDER!
      pd.basnoff = [pr.basnoff pr.N+pr.t.basnoff]; % effective setupbasisdofs
      pd.N = pr.N + pr.t.N;            % finish up appending to pd basis setup

      % Evaluate all bases, 3 styles (no periodic direct image sums here)...
      if nargout==1, A = pd.evalbases@problem(p, opts);
      elseif nargout==2, [A Ax] = pd.evalbases@problem(p, opts);
      else [A Ax Ay] = pd.evalbases@problem(p, opts); end
      
      % Set up for image sums of obst only, only in the untranslated extdom...
      extdomi = find([pr.doms.isair]);          % index of the extdom
      if numel(extdomi)~=1, error('there appears to be >1 exterior domain!');end
      i = pr.domainindices(p)==extdomi;         % logical array, size of p.x
      opts.dom = pr.extdom;                     % speeds up prob.evalbases loop
      j = 1:pr.N;                               % obst dof indices
      pc = pointset(p.x(i));                    % only eval images where need!
      if ~isempty(p.nx), pc.nx = p.nx(i); end   % need matching normals too
      for n=[-pr.nei:-1 1:pr.nei]               % direct sum nei contribs
        pc.x = p.x(i) - pr.d*n; ph = pr.a^n;    % target translation, phase
        if nargout==1, A(i,j) = A(i,j) + ph * pr.evalbases@problem(pc, opts);
        elseif nargout==2, [Ac Axc] = pr.evalbases@problem(pc, opts);
          A(i,j) = A(i,j) + ph * Ac; Ax(i,j) = Ax(i,j) + ph * Axc;
        else [Ac Axc Ayc] = pr.evalbases@problem(pc, opts);
          A(i,j) = A(i,j) + ph * Ac; Ax(i,j) = Ax(i,j) + ph * Axc;
          Ay(i,j) = Ay(i,j) + ph * Ayc; end
      end
    end

    function showfullfield(pr, o) % ................... overloads @scattering
    % SHOWFULLFIELD - eval and plot u_i+u on grid (Re part), periodic wrapping
    %
    % By default, shows width of 3 periods, and wraps so u looks quasiperiodic
    %
    %   opts.imag = true, plots imag instead of real part
    %   opts.bdry = true, shows boundary too (& options passed to showbdry)
    %   opts.nowrap = true, goes back to default bounding-box (no wrapping).
    %   opts.nx = # gridpoints across one period
    %   opts.ymax = sets y-range to [-ymax ymax]
   
      if nargin<2, o = []; end
      if ~isfield(o, 'imag'), o.imag = 0; end
      if ~isfield(o, 'bdry'), o.bdry = 0; end
      if ~isfield(o, 'nowrap'), o.nowrap = 0; end
      if ~isfield(o, 'nx'), o.nx = 50; end
      if ~isfield(o, 'ymax'), o.ymax = 0; end

      if ~o.nowrap                   % bb make one period, then wrap + expand
        if o.ymax==0       % get bb from obstacle geom in default way
          oo = pr.gridboundingbox; o.bb = oo.bb; clear oo;
        else, o.bb(3) = -o.ymax; o.bb(4) = o.ymax; end
        o.dx = pr.d/o.nx; o.bb(1) = -pr.d/2+o.dx/2; % overwrite dx and x-range
        tiny = 1e-12; o.bb(2) = pr.d/2-o.dx/2+tiny;
        u = pr.gridsolution(o);                     % eval scatt field (slow)
        [ui gx gy di] = pr.gridincidentwave(o); u = ui + u; % total field
        u = [pr.a^(-1)*u u pr.a*u];   % three phased copies of total field
        di = [di di di];              % in case needed (isn't as of yet)
        o.bb(1) = o.bb(1) - pr.d; o.bb(2) = o.bb(2) + pr.d; % widen the bb
        gx = [gx-pr.d gx gx+pr.d];              % needed for plotting
      else, u=pr.gridsolution(o); % let bb be chosen by default elsewhere
        [ui gx gy di] = pr.gridincidentwave(o); u=ui+u;
      end % compute tot field
        
      figure;
      if o.imag, imagesc(gx, gy, imag(u));title('Im[u_{tot}]');
      else, imagesc(gx, gy, real(u)); title('Re[u_{tot}]');
      end
      utils.goodcaxis(u);
      axis equal tight; colorbar; set(gca,'ydir','normal'); hold on;
      if o.bdry, pr.showbdry(struct('nei',1)); end
    end      
 
    function h = showbdry(pr, o)   % ........................ crude plot bdry
    % SHOWBDRY - shows boundary segments in a qpscatt, with their natural sense
    %
    % h = SHOWBDRY(pr, opts) passed in various options, including
    %   opts.nei = 0,1,... how many periodic copies either side to plot
      if nargin<2, o = []; end
      if ~isfield(o, 'nei'), o.nei = 1; end       % default # neighbors to plot
      h = []; on = ones(size(pr.segs));
      for n=-o.nei:o.nei
        h = [h; domain.showsegments(pr.segs.translate(n*pr.d), on, o)];
      end
    end
    
  end % methods
end
