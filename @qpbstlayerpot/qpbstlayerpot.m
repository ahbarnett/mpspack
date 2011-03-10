% QPBSTLAYERPOT - layer potential basis & phased copies for QP strip unit cell
%
% b = QPBSTLAYERPOT(st, s, clist, a) creates set of copies of layer densities
%  on segment s, 
%
% cut down from qpuclayerpot
%
% See also: BOUNDEDQPSTRIP

classdef qpbstlayerpot < handle & basis
  properties
    st                     % the bounded qp strip domain (needed for st.e vec)
    a                      % 1-by-2, mixture weights of SLP and DLP
    lp                     % layer-potential object on principal seg
    clist                  % list of copy translations in multiples of st.e vec
    quad                   % quadrature option (opts.quad in S,D,T)
  end
  
  methods
    function b = qpbstlayerpot(st, s, clist, a, opts) % ............ constructor
      if nargin<5, opts = []; end
      b.st = st;
      if ~isnumeric(a)
        switch a
         case {'single', 'S', 's', 'SLP'}
          a = [1 0];
         case {'double', 'D', 'd', 'DLP'}
          a = [0 1];
        end
      end
      b.a = a;
      b.lp = layerpot(s, a, opts);
      b.clist = clist;
      b.lp.doms = st;   % make the LP affect the strip domain
    end % func
    
    function Nf = Nf(b, opts)  %...............Nf method reads segment # quadr
      Nf = numel(b.lp.seg.x);
    end
    
    function showgeom(b, varargin) % ............. crude show discr pts of segs
      for j=1:numel(b.clist)
        seg(j) = b.lp.seg.translate(b.clist(j)*b.st.e); % temp copies to plot
      pm(j) = 1; end                             % dummy pm choices
      domain.showsegments(seg, pm, varargin{:}); % pass opts if present onwards
    end
    
    function [A A1 A2] = eval(b, p, opts); % ......... evaluator matrix
      % really this stuff should be universal to all basis sets, as evalcopies
      if nargin<3, opts = []; end
      a = b.st.a;                         % alpha Bloch phase
      A = zeros(numel(p.x), numel(b.lp.seg.x)); A1 = A; A2 = A;
      ph = a.^b.clist;                  % phases
      for j=1:numel(b.clist)
        if b.clist(j)==0, pj = p;      %  untranslated: preserve identity of seg
        else
          pj = pointset(p.x - b.clist(j)*b.st.e, p.nx); % or transl target pts
        end
        if nargout==1
          A = A + ph(j) * b.lp.eval(pj, opts);
        elseif nargout==2
          [Aj A1j] = b.lp.eval(pj, opts);      % phased copies
          A = A + ph(j)*Aj; A1 = A1 + ph(j)*A1j;
        elseif nargout==3
          [Aj A1j A2j] = b.lp.eval(pj, opts);
          A = A + ph(j)*Aj; A1 = A1 + ph(j)*A1j; A2 = A2 + ph(j)*A2j;
        end
      end
    end
    
    function Q = evalboundedstripdiscrep(b, s, opts)
    % EVALBOUNDEDSTRIPDISCREP - matrix mapping basis coeffs to discrepancy
    %
    % Special routine for the QPBSTLAYERPOT basis object which is tied to a
    % bounded qp strip unit-cell; omits terms which cancel (buf=0 always).
    %
    % Currently for L&R, case of clist=[0 1] and [-1 2] are hacked - fix?
    % Also codes explicit cancellation in C-type blocks for general LPs.
      if ~isa(s, 'boundedqpstrip')
        error('s must be a boundedqpstrip object!'); end
      if nargin<3, opts = []; end
      d = s.e;                               % strip width (s taken over b.st)
      
      if s.buf==0 & numel(b.clist)==2 & b.clist==[0 1] & b.lp.seg==s.L
        N = numel(s.L.x);  % hack L&R QP 1-box scheme, expl jump term...
        pr = pointset(s.L.x + b.st.e, s.L.nx); % R pts
        pll = pointset(s.L.x - b.st.e, s.L.nx); % copy of pts one left of L
        [Q Qn] = b.lp.eval(pr, opts);
        [Ql Qnl] = b.lp.eval(pll, opts);
        a = s.a;           % Bloch phase from qpstrip domain
        Q = [b.a(2)*eye(N); -b.a(1)*eye(N)] + [a*Ql-Q/a; a*Qnl-Qn/a];
        
      elseif s.buf==1 & numel(b.clist)==2 & b.clist==[-1 2] & b.lp.seg==s.L
        fprintf('qpbstlayerpot.evalboundedstripdiscrep: using 3-box scheme\n');
        N = numel(s.L.x);  % hack L&R QP 3-box scheme, expl jump term...
        pr = pointset(s.L.x + 3*b.st.e, s.L.nx); % R pts
        pll = pointset(s.L.x - 3*b.st.e, s.L.nx); % copy of pts one left of L
        [Q Qn] = b.lp.eval(pr, opts);
        [Ql Qnl] = b.lp.eval(pll, opts);
        a = s.a^3;           % Bloch phase across 3x qpstrip domain
        Q = [b.a(2)*eye(N); -b.a(1)*eye(N)] + [a*Ql-Q/a; a*Qnl-Qn/a];
        
      elseif s.buf==0 & numel(b.clist)==b.clist(end)-b.clist(1)+1 % contig list
        % hack nearby cancellation in C-type blocks explicitly ...
        nr = b.clist(1)-1; nl = b.clist(end); % translations, also alpha-powers
        pr = pointset(s.L.x - nr*b.st.e, s.L.nx); % shifted R target 
        pl = pointset(s.L.x - nl*b.st.e, s.L.nx); % shifted L target
        [Q Qn] = b.lp.eval(pr, opts);
        [Ql Qnl] = b.lp.eval(pl, opts);
        a = s.a;           % Bloch phase from qpstrip domain
        Q = [Ql*a^nl-Q*a^nr; Qnl*a^nl-Qn*a^nr];  % alpha-powers
        
      else           % general discrep computation... cancellation is numerical
        Q = evalboundedstripdiscrep@basis(b, s, opts);
      end
    end
    
  end % methods
end
