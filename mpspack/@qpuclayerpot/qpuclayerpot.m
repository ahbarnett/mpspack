classdef qpuclayerpot < handle & basis

% QPUCLAYERPOT - create sticking-out copies of layer densities for a unit cell
%
%  This creates 
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
    Pmp Pop Ppp Opo Opp Pmm Pom Ppm Omo Omp  % Bloch-k-indep storage matrices
                                             % for Q recalcs, key: p=+1 o=0 m=-1
    k_storage                       % k for which stored P??, O?? were computed
  end
  
  methods
    function b = qpuclayerpot(uc, direc, a, k, opts)
    % QPUCLAYERPOT - setup a sticking-out QP layer potential basis on unit cell
    %
    %  bas = qpuclayerpot(uc, direc, a, k, opts) creates a set of six
    %   copies of density on either L or B seg of uc, depending on if direc is
    %   'L' or 'B'. a controls mix of SLP and DLP as in layerpot.
    %  Note that in fact the LPs are not constructed until evaluation time, but
    %   the setup of locations and phases is done here.
    %
    % See also: LAYERPOT, QPUNITCELL
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
    %  automatically knows which L or B orientation
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
    
    function Q = evalunitcelldiscrep(b, uc, opts) % ...... overloads from basis
    % EVALUNITCELLDISCREP - return Q matrix mapping a QPUC basis to UC discrep
    %
    %  Does include source quadr wei, ie dofs are (density value_j)
    %  Currently uses no symmetry relations. Stores Bloch-indep interaction
    %   matrices in b.P??, b.O??.
    %  Uses only distant interactions + Id, so has spectral convergence.
    %  See testqpuclayerpot.m for convergence plot
    %
    % Also see: QPUCLAYERPOT, LAYERPOT, QPUNITCELL
    
      % static storage for interaction rep matrices... needs know about uc.bas
      b.updateNf;
      k = b.k; a = b.a; p = b.pseg; o = b.oseg; % p,o: parallel & other segs
      if isempty(b.Pmp) | b.k~=b.k_storage | b.Nf~=size(b.Pmp, 2) % recompute?
        b.k_storage = b.k;         % remember for what k P?? O?? will be for
        % compute matrices by creating source seg copy c & moving it around
        % (note, utils.copy in new segment create is slow, only do once...)
        c = p.translate(-b.ep+b.eo); l = layerpot(c, a, k); % parallel seg copy
        [b.Pmp t] = l.eval(p); b.Pmp = [b.Pmp; t]; % stack n-derivs below vals
        c.translate(b.ep); l = layerpot(c, a, k); % c = b.eo
        [b.Pop t] = l.eval(p); b.Pop = [b.Pop; t];
        c.translate(b.ep); l = layerpot(c, a, k); % c = b.ep + b.eo
        [b.Ppp t] = l.eval(p); b.Ppp = [b.Ppp; t];
        c.translate(-b.eo); l = layerpot(c, a, k); % c = b.ep
        [b.Opo t] = l.eval(o); b.Opo = [b.Opo; t];
        c.translate(b.eo); l = layerpot(c, a, k);  % c = b.ep + b.eo
        [b.Opp t] = l.eval(o); b.Opp = [b.Opp; t];
        % these should probably be replaced by transposes / reflections ...
        c.translate(-2*(b.ep+b.eo)); l = layerpot(c, a, k); % c = -b.ep - b.eo
        [b.Pmm t] = l.eval(p); b.Pmm = [b.Pmm; t];
        c.translate(b.ep); l = layerpot(c, a, k); % c = -b.eo
        [b.Pom t] = l.eval(p); b.Pom = [b.Pom; t];
        c.translate(b.ep); l = layerpot(c, a, k); % c = b.ep -b.eo
        [b.Ppm t] = l.eval(p); b.Ppm = [b.Ppm; t];
        c.translate(-3*b.ep + b.eo); l = layerpot(c, a, k); % c = -2*b.ep
        [b.Omo t] = l.eval(o); b.Omo = [b.Omo; t];
        c.translate(b.eo); l = layerpot(c, a, k); % c = -2*b.ep + b.eo
        [b.Omp t] = l.eval(o); b.Omp = [b.Omp; t];
      end
      p = b.facp(b.uc); o = b.faco(b.uc);       % parallel & other, phase facs
      %fprintf('p,o = %g,%g\n', p,o)
      % Note we're stacking f, f' together here, stored hit by phase factors...
      Qp = (o/p)*b.Pmp + o*b.Pop + o*p*b.Ppp - (1/o/p)*b.Pmm - (1/o)*b.Pom - (p/o)*b.Ppm;
      Qo = p*b.Opo + p*o*b.Opp - (1/p^2)*b.Omo - (o/p^2)*b.Omp;
      Qp = Qp + [a(2)*eye(b.Nf); -a(1)*eye(b.Nf)];  % jump relations (not 1/2)
      if b.direc=='L'          % make sure f & g are correct way round
        Q = [Qp; Qo];          % f,f' hit by parallel contribs
      else
        Q = [Qo; Qp];          % g,g' hit by parallel contribs
      end
    end % func
  end % methods
end
