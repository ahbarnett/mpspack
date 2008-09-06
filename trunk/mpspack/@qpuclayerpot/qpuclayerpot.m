    % QPUCLAYERPOT - setup a sticking-out QP layer potential basis on unit cell
    %
    %  bas = qpuclayerpot(uc, direc, a, k) creates a set of six copies of a
    %   layer density on either L or B segment of the unit cell object uc,
    %   depending on if direc is
    %   'L' or 'B'. The argument a controls mix of SLP and DLP, as in layerpot.
    %  Note that in fact the LPs are not constructed until evaluation time, but
    %   the setup of locations and phases is done here. Similarly, the phases
    %   used for the copies depend on (alpha, beta) for the unit cell at the
    %   time of evaluation.
    %
    % See also: LAYERPOT, QPUNITCELL

classdef qpuclayerpot < handle & basis
  properties
    uc                              % unit cell that it's attached to
    pseg                            % original source segment (will be L or B)
    oseg                            % other (non-parallel) source segment (B/L)
    direc                           % 'L' or 'B'
    ep, eo                          % lattice vectors in parallel & other direc
    facp, faco                      % phase factors in parallel & other direc
    segf                            % eval segment list func (first is L or B)
    phasef                          % eval phase factor func @(parallel, other)
    a                               % 1-by-2, mixture weights of SLP and DLP
    quad                            % quadrature option (opts.quad in S,D,T)
  end
  
  methods
    function b = qpuclayerpot(uc, direc, a, k) % ................ constructor
      b.uc = uc;
      if direc=='L' | direc=='l'
        s = uc.L; b.oseg = uc.B; b.direc = 'L';
        b.ep = uc.e2; b.eo = uc.e1; p = @(uc) uc.b; o = @(uc) uc.a; 
      elseif direc=='B' | direc=='b'
        s = uc.B; b.oseg = uc.L; b.direc = 'B';
        b.ep = uc.e1; b.eo = uc.e2; p = @(uc) uc.a; o = @(uc) uc.b;
      else
        error('invalid direc!');
      end
      b.pseg = s; b.facp = p; b.faco = o;
      b.updateNf;
      % Note that the following segs and phases are for direct eval only...
      b.segf = @(ep,eo) [s translate(s, ep) translate(s, -ep) ...
                         translate(s, eo) translate(s, eo+ep) ...
                         translate(s, eo-ep)]; % seg list func, 1+5 copies
      b.phasef = @(uc) [1 p(uc) 1./p(uc) o(uc) o(uc).*p(uc) o(uc)./p(uc)];
      if ~isnumeric(a)
        switch a
         case {'single', 'S', 's', 'SLP'}
          a = [1 0];
         case {'double', 'D', 'd', 'DLP'}
          a = [0 1];
        end
      end
      b.a = a; b.k = k;
    end % func
    
    function updateNf(b)  % ............... Nf property reads segment # quadr
      b.Nf = numel(b.pseg.x);      % # quadr points on 'parallel' segment
    end
    
    function showgeom(b, opts) % .................. crude show discr pts of seg
      seg = b.segf(b.ep, b.eo);          % make temporary copy of segs for plot
      domain.showsegments(seg);
    end
    
    function [A Ax Ay] = eval(b, p, opts) % ...... standard distant evaluation
    % EVAL - standard evaluation of QP unit-cell layer potential copies basis
    %
    %  Refers to current alpha, beta and segments in unit cell, but nice that
    %  automatically knows which L or B orientation it was set up with.
      b.updateNf;
      seg = b.segf(b.ep, b.eo);                   % construct seg copies
      phase = b.phasef(b.uc);
      A = zeros(numel(p.x), b.Nf);
      if nargout>=1, Ax = A; end
      if nargout>=2, Ay = A; end
      for j=1:numel(seg)
        lp = layerpot(seg(j), b.a, b.k);
        ph = phase(j);
        if nargout==1
          A = A + ph * lp.eval(p);      % could add jump rels if p=seg
        elseif nargout==2
          [As Asx] = lp.eval(p);
          A = A + ph * As; Ax = Ax + ph * Asx;
        else
          [As Asx Asy] = lp.eval(p);
          A = A + ph * As; Ax = Ax + ph * Asx; Ay = Ay + ph * Asy;
        end
      end
    end
    
    function [varargout] = evaltargetcopies(b, p, uc, opts) % ...overloads basis
    % EVALTARGETCOPIES - see basis/EVALTARGETCOPIES for interface
    %
    %  For this qpuclayerpot basis type, a special copy list is sent to generic
    %   basis routine. Ie, it's a wrapper around the layerpot routine.
    %   Note that target-movement is used to achieve source move for nei=0
      if nargin<4, opts = []; end
      if isfield(opts, 'nei'), nei = opts.nei; else; nei = 0; end
      if nei~=0 & nei~=1, error('opts.nei must be 0 or 1!'); end
      if nei==0              % special copylist, six LP sources with phases...
        x = uc.e1; y = uc.e2; be = uc.b;
        if b.direc=='L'      % 'L' direction (left)
          c.t = [y, 0, -y, -x+y, -x, -x-y];  % -ve since faking w/ targ move
          c.apow = [0, 0, 0, 1, 1, 1];
          c.remph = [be^-1, 1, be, be^-1, 1, be];
        else                 % 'B' direction (bottom)
          c.t = [x, 0, -x, -y+x, -y, -y-x]; % x<->y swapped
          c.apow = [-1, 0, 1, -1, 0, 1];
          c.remph = [1, 1, 1, be, be, be];
        end
        opts.copylist = c;
        lp = layerpot(b.pseg, b.a, b.k, opts);
        [varargout{1:nargout}] = lp.evaltargetcopies@basis(p, uc, opts);
      else                   % nei=1, in which case special copies list needed
        x = uc.e1; y = uc.e2; be = uc.b;
        if b.direc=='L'      % 'L' direction (left)
          c.t = [x+2*y, x, x-2*y, -2*x+2*y, -2*x, -2*x-2*y];
          c.apow = [-1, -1, -1, 2, 2, 2];
          c.remph = [be^-2, 1, be^2, be^-2, 1, be^2];
        else                 % 'B' direction (bottom)
          c.t = [y+2*x, y, y-2*x, -2*y+2*x, -2*y, -2*y-2*x]; % x<->y swapped
          c.apow = [-2, -0, 2, -2, 0, 2];
          c.remph = [be^-1, be^-1, be^-1, be^2, be^2, be^2];          
        end
        opts.copylist = c;
        lp = layerpot(b.pseg, b.a, b.k, opts);
        % NB to test `naive' sum, change lp to b here, & kill opts.copylist...
        [varargout{1:nargout}] = lp.evaltargetcopies@basis(p, uc, opts);
      end
    end
    
    function [Q d] = evalunitcelldiscrep(b, uc, opts) % ....overloads from basis
    % EVALUNITCELLDISCREP - return Q matrix mapping a QPUC basis to UC discrep
    %
    %  Same usage as basis/EVALUNITCELLDISCREP, which it overloads.
    %  (Includes source quadr wei, ie dofs are density values, as always).
    %  Currently uses no symmetry relations. Stores Bloch-indep interaction
    %   matrices in d.P??, d.O?? and their k in d.k (changed d from b 8/20/08)
    %  Uses only distant interactions + Id, so has spectral convergence.
    %  See testqpuclayerpot.m for convergence plot
    %
    % See also QPUCLAYERPOT, LAYERPOT, QPUNITCELL, basis/EVALUNITCELLDISCREP.
      uc = b.uc;                            % unit cell in bas takes priority
      if nargin<3, opts = []; end
      if isfield(opts, 'nei') & opts.nei~=0, error('opts.nei must be 0!'); end
      if isfield(opts, 'data'), d = opts.data; else d = []; end
      k = uc.k;                                   % wavenumber (UC beats bas)
      b.updateNf;
      recompute = isempty(d);
      if ~recompute              % make sensible guess as to if needs bas evals
        recompute = (k~=d.k | b.Nf~=size(d.Pmp,2));
      end
      wantpoly = 0; if isfield(opts, 'poly') & opts.poly, wantpoly = 1; end
      a = b.a;                   % SLP, DLP mixture coeffs
      if recompute
        d.k = k;         % remember for what wavenumber data P?? O?? is for
        p = b.pseg; o = b.oseg; % p,o: 'parallel' & 'other' segments
        % compute matrices by creating source seg copy c & moving it around
        % (note, utils.copy in new segment create is slow, only do once...)
        c = p.translate(-b.ep+b.eo); l = layerpot(c, a, k); % parallel seg copy
        [d.Pmp t] = l.eval(p); d.Pmp = [d.Pmp; t]; % stack n-derivs below vals
        c.translate(b.ep); l = layerpot(c, a, k); % c = b.eo
        [d.Pop t] = l.eval(p); d.Pop = [d.Pop; t];
        c.translate(b.ep); l = layerpot(c, a, k); % c = b.ep + b.eo
        [d.Ppp t] = l.eval(p); d.Ppp = [d.Ppp; t];
        c.translate(-b.eo); l = layerpot(c, a, k); % c = b.ep
        [d.Opo t] = l.eval(o); d.Opo = [d.Opo; t];
        c.translate(b.eo); l = layerpot(c, a, k);  % c = b.ep + b.eo
        [d.Opp t] = l.eval(o); d.Opp = [d.Opp; t];
        % these should probably be replaced by transposes / reflections ...
        c.translate(-2*(b.ep+b.eo)); l = layerpot(c, a, k); % c = -b.ep - b.eo
        [d.Pmm t] = l.eval(p); d.Pmm = [d.Pmm; t];
        c.translate(b.ep); l = layerpot(c, a, k); % c = -b.eo
        [d.Pom t] = l.eval(p); d.Pom = [d.Pom; t];
        c.translate(b.ep); l = layerpot(c, a, k); % c = b.ep -b.eo
        [d.Ppm t] = l.eval(p); d.Ppm = [d.Ppm; t];
        c.translate(-3*b.ep + b.eo); l = layerpot(c, a, k); % c = -2*b.ep
        [d.Omo t] = l.eval(o); d.Omo = [d.Omo; t];
        c.translate(b.eo); l = layerpot(c, a, k); % c = -2*b.ep + b.eo
        [d.Omp t] = l.eval(o); d.Omp = [d.Omp; t];
      end
      if wantpoly % ....... separate matrix contribs by alpha powers
        ML = numel(uc.L.x); MB = numel(uc.B.x); % if L,R diff # pts from B,T?
        Q = zeros(2*(ML+MB), b.Nf, 5);
        be = uc.b;                     % beta
        if b.direc=='L'                % al is 'other', be is 'parallel'
          Q(1:2*ML,:,2) = -(1/be)*d.Pmm - d.Pom - be*d.Ppm;            % al^{-1}
          Q(:,:,3) = [a(2)*eye(b.Nf); -a(1)*eye(b.Nf); ...
                      be*d.Opo - (1/be^2)*d.Omo];                      % al^0
          Q(:,:,4) = [(1/be)*d.Pmp + d.Pop + be*d.Ppp; ...
                      be*d.Opp - (1/be^2)*d.Omp];                      % al^1
        else                           % al is 'parallel', be is 'other'
          Q(1:2*ML,:,1) =  -d.Omo - be*d.Omp;                          % al^{-2}
          Q(2*ML+1:end,:,2) =  be*d.Pmp -(1/be)*d.Pmm;                 % al^{-1}
          Q(2*ML+1:end,:,3) =  be*d.Pop -(1/be)*d.Pom + ...
              [a(2)*eye(b.Nf); -a(1)*eye(b.Nf)];                       % al^0
          Q(:,:,4) = [d.Opo + be*d.Opp; be*d.Ppp - (1/be)*d.Ppm];      % al^1
        end
      else     % ....... compute all al, be powers to get single Q matrix
        % look up parallel & other, phase facs... (MUCH faster than b.facp(uc))
        if b.direc=='L', p = uc.b; o = uc.a; else p = uc.a; o = uc.b; end 
        %fprintf('p,o = %g,%g\n', p,o)
        % Note we're stacking f, f' together here, stored hit by phase factors
        Qp = (o/p)*d.Pmp + o*d.Pop + o*p*d.Ppp - (1/o/p)*d.Pmm - (1/o)*d.Pom - (p/o)*d.Ppm;
        Qo = p*d.Opo + p*o*d.Opp - (1/p^2)*d.Omo - (o/p^2)*d.Omp;
        Qp = Qp + [a(2)*eye(b.Nf); -a(1)*eye(b.Nf)];  % jump relations (not 1/2)
        if b.direc=='L'          % make sure f & g are correct way round
          Q = [Qp; Qo];          % f,f' hit by parallel contribs
        else
          Q = [Qo; Qp];          % g,g' hit by parallel contribs
        end
      end
    end % func
  end % methods
end
