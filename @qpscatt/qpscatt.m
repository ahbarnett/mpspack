% QPSCATT - define quasi-periodic (grating) frequency-domain scattering problem
%
%  pr = qpscatt(airdoms, doms, d) creates a scattering problem object pr
%   as with scattering, except a quasi-periodic version with x-periodicity of d.
%   The solution method will be via FTyLPs.
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

% (C) 2010 Alex Barnett

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
    kpyfix                                  % list of k_y for poles to cross
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
        setincidentwave@scattering(pr, t); % call superclass method
        om = pr.k; kvec = om*exp(1i*t); pr.a = exp(1i*real(conj(kvec) * pr.d));
        pr.t.setbloch(exp(1i*real(conj(kvec) * pr.t.e))); % set Bloch alpha
        % set the k-poles list in pr object...
        n = 2*ceil(om*pr.d/2/pi)+10;  % # diffraction orders?
        kpx = om*cos(t)+(-n:n)*2*pi/pr.d; kpy = sqrt(om^2-kpx.^2); % kp y-compt
        [dum i] = sort(abs(kpy),'ascend');       % note +ve sqrt key here
        pr.kpx = kpx(i); pr.kpy = kpy(i);
        % reset the QP FTyLP bases, and add any Wood's anomaly fix bases...
        for i=1:2, pr.t.bas{i}.requadrature(pr.t.bas{i}.N, ...
                      struct('omega',om, 'nearsing', abs(pr.kpy(1)))); end
        pr.t.bas = pr.t.bas(1:2); pr.kpyfix = []; % kill off any Wood fix bases
      end
    end
    
    function showkyplane(pr)
    % SHOWKYPLANE - plot Sommerfeld quadrature nodes, k_y poles, +-om, etc
      figure; plot(pr.t.bas{1}.lp.kj, '+'); axis tight equal; hold on;
      plot([pr.kpy -pr.kpy], 'rx');
      plot(real(pr.kpyfix),imag(pr.kpyfix),'go');
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
      rhs = [obstrhs; zeros(pr.t.N,1)];    % QP block square, so use # QP dofs
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
      Mq = pr.t.bas{1}.N;              % # FTyLP quadr pts
      nb = numel(pr.t.bas);            % how many QP bases in t
      nei = pr.nei; buf = pr.buf; d = pr.d; a = pr.a; % get Bloch alpha
      
      % fill Q...
      Q = t.evalbasesdiscrep();        % is block-diag - need not really fill!
      
      % fill A: self, then directly sum neighbor images contribs...
      A = fillbcmatrix@problem(pr);    % obstdofs->mismatch block
      %for i=1:numel(bob.bas)                 % all bases in problem or basis obj
       % b = bob.bas{i}; ns = bob.basnoff(i)+(1:b.Nf); % dof indices for this bas  
      %for n=[-nei:-1 1:nei], A = A + a^n * l.eval(pointset(s.x-d*n)); end

      % TO FINISH!!! - use fillbcmatrix@problem with transl and o.doms
      % restrictions!
      
      
      
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
      C = zeros(Nq, N); f = utils.copy(t.bas{1}.lp); %dummy bas for transl evals
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
      
      pr.A = [A B; C Q];      % stack the E matrix for output (it's called A)
            
      if ~isempty(pr.kpyfix), % augment matrix with Wood fix ...
      end
      
    end % fillbcmatrix
    
    function [u di] = pointsolution(pr, p) % ............... overloads @problem
    % POINTSOLUTION - evaluate solution to a problem on a pointset, given coeffs
    %
    % Overloads problem.pointsolution, but creating dumm problem and dummy
    %  extdom which marries the bases from the obst extdom and the qpstrip
    %
    % See also PROBLEM.POINTSOLUTION
    
      % make dummy problem with new extdom object (whose d.bas are as before)
      pd = utils.copy(pr); pd.A = []; % (a hack!). leave p.co since need it
      for i=1:numel(pd.doms), d = pr.doms(i); if d.isair, 
          d = utils.copy(d); pd.doms(i) = d; pd.extdom = d; end, end
      % now add the QP bases onto those in obst exterior domain...
      pd.extdom.bas = {pd.extdom.bas{:} pd.t.bas{:}}; % marry the bases lists
      pd.setupbasisdofs;
      [u di] = pointsolution@scattering(pd, p);  % pass to usual evaluator
    end

   function showfullfield(pr, o) % ................... overloads @scattering
    % SHOWFULLFIELD - eval and plot u_i+u on grid (Re part), periodic wrapping
    %
    % By default, shows width of 3 periods, and wraps so u quasiperiodic
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
        o.dx = pr.d/o.nx; o.bb(1) = pr.t.Lo+o.dx/2; % overwrite dx and x-range
        tiny = 1e-12; o.bb(2) = pr.t.Ro-o.dx/2+tiny;
        u = pr.gridsolution(o);                     % eval scatt field (slow)
        [ui gx gy di] = pr.gridincidentwave(o); u = ui + u; % total field
        u = [pr.a^(-1)*u u pr.a*u];   % three phased copies of total field
        di = [di di di];              % in case needed (isn't as of yet)
        o.bb(1) = o.bb(1) - pr.d; o.bb(2) = o.bb(2) + pr.d; % widen the bb
        gx = [gx-pr.d gx gx+pr.d];              % needed for plotting
      else, u=pr.gridsolution(o); % let bb be chosen by default elsewhere
        [ui gx gy di] = pr.gridincidentwave(o); u=ui+u; end % compute tot field
        
      figure;
      if o.imag, imagesc(gx, gy, imag(u));title('Im[u_{tot}]');
      else, imagesc(gx, gy, real(u)); title('Re[u_{tot}]');
      end
      c = caxis; caxis([-1 1]*0.7*max(c));
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
