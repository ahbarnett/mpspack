% BOUNDEDQPSTRIP - create a quasi-periodic strip domain of bounded y extent
%
% st = boundedqpstrip([B T], k) creates a strip bounded by segments B at the
%   bottom (oriented in +x sense), and T at the top (-x sense). B and T must
%   have the same x extents, so that left and right walls can be vertical.
%   k is the wavenumber in the strip.
%
% st = qpstrip([B T], k, opts) includes the following options:
%    opts.M = sets # quadr nodes on left & right walls
%    opts.quad = choose quadrature type for L & R walls
%    opts.buf = buffer width of strip (0 gives usual, 1 gives 3xwidth, etc)
%
% ISSUES:
%  * changed so B seg can be reversed, and the creator will deal with it.
%
% See also DOMAIN, QPSTRIP

% (C) 2010 Alex Barnett
classdef boundedqpstrip < handle & domain
  properties
    a                             % alpha QP phase factor
    e                             % horizontal lattice vector (complex #)
    r                             % reciprocal lattice vector (complex #)
    Lo, Ro                        % left & right x-limits
    B, T, L, R                    % segments: B,T, L, R walls
    buf                           % buffer width of strip (0,1,...)
    N                    % total # degrees of freedom in QP basis sets
    basnoff              % dof offsets for QP basis sets (numeric list)
  end
  
  methods
    function st = boundedqpstrip(s, k, o) % ....................... constructor
      if nargin<3, o = []; end
      if ~isfield(o, 'M'), o.M = 20; end
      if ~isfield(o, 'quad'), o.quad = 'g'; end
      if ~isfield(o, 'buf'), o.buf = 0; end
      if ~isa(s, 'segment')
        error('1st argument must be array of segment objects!'); end
      B = s(1); T = s(2);
      Lo = real(B.eloc(1)); Ro = real(B.eloc(2)); pmB = 1; % default pm for B
      if Lo>Ro, temp = Ro; Ro = Lo; Lo = temp; pmB = -1;
        warning('swapped sense of B seg...');
      end
      e = Ro-Lo;
      L = segment(o.M, [B.eloc(1+(1-pmB)/2) T.eloc(2)], o.quad);
      R = translate(L, e);
      tiny = 1e-14;
      if abs(real(T.eloc(1))-Ro)+abs(real(T.eloc(2))-Lo)>tiny
        warning('segments do not appear to match in x-limits!');
      end
      st = st@domain([L B R T], [-1 pmB 1 1]);     % use bounded domain (as uc)
      st.L = L; st.R = R; st.e = e; st.B = B; st.T = T; st.Lo = Lo; st.Ro = Ro;
      st.setbloch;
      st.k = k;
      st.r = 2*pi/st.e;
      st.buf = o.buf;
    end

    function setbloch(st, a)
    % SETBLOCH - sets alpha phase for QP strip object
      if nargin<2, a = 1; end    % default
      st.a = a;
    end
    
    function [N noff] = setupbasisdofs(uc) % .................................
    % SETUPBASISDOFS - set up indices of basis degrees of freedom in unit cell
    %
    % copied from qpunitcell.
    if isempty(uc.bas), fprintf('warning: no basis sets in unit cell!\n'); end
      noff = zeros(1, numel(uc.bas));
      N = 0;                         % N will be total # dofs (cols of B)
      for i=1:numel(uc.bas)
        noff(i) = N; N = N + uc.bas{i}.Nf; end    % set up dof offsets of bases
      uc.N = N; uc.basnoff = noff;   % store stuff as unit cell properties
    end

    function addqpftylayerpots(st, o) % ...............................
    % ADDQPLAYERPOTS - add SLP+DLP layer potentials to L in a boundedqpstrip
    %
    % All calling opts are passed to layerpot, apart from opts.omega which
    % is inherited from the boundedqpstrip domain.
      if nargin<2, o = []; end
      if st.buf>0, error 'buf>0 not implemented with FTyLPs'; end
      if ~isfield(o,'omega'), o.omega = st.k; end      % inherit wavenumber
      % compute distance from origin to nearest singularity given Bloch params:
      o.nearsing=min(abs(sqrt(st.k^2-(log(st.a)/1i+(-100:100)*2*pi/st.e).^2)));
      st.bas = {st.bas{:}, ftylayerpot(st.Lo, 'd', o), ...
                ftylayerpot(st.Lo, 's', o)};
      for i=0:1, st.bas{end-i}.doms = st; end  % make strip the affected domain
      st.setupbasisdofs;
    end
    
    function addqpbstlayerpots(t, o) % ..................................
    % ADDQPBSTLAYERPOTS - add left-right QP scheme LPs to bounded strip domain
    %
    % Note: order was swapped to DLP then SLP, to keep diagonal in right place
    % Barnett 2/15/11
      if nargin<2, o = []; end
      clist = [-t.buf, 1+t.buf]; % copylist for source segments
      t.bas{1} = qpbstlayerpot(t, t.L, clist, 'd', o); t.bas{1}.doms = t;
      t.bas{2} = qpbstlayerpot(t, t.L, clist, 's', o); t.bas{2}.doms = t;
      t.setupbasisdofs;
    end
    
    function requadrature(uc, N) % ....... overloads domain.requadrature
    % REQUADRATURE - reset # quadrature pts for L, R wall layerpots
      if nargin<2, N = []; end       % uses default
      uc.L.requadrature(N); uc.R.requadrature(N);
      uc.setupbasisdofs;
    end
    
    function Q = evalbasesdiscrep(st, opts) % .................... Q
    % EVALBASESDISCREP - fill Q matrix mapping strip bases coeffs to discrepancy
    %
    %  Q = EVALBASESDISCREP(st) returns matrix of basis function evaluations
    %   given basis set dofs as stored in the bounded strip domain.
    %
    % Based on: qpunitcell.evalbasesdiscrep
      if nargin<2, opts = []; end
      N = st.N; noff = st.basnoff;    % # basis dofs
      M = 2*numel(st.L.x);            % # discrep dofs (f and f')
      Q = zeros(M,N);                 % preallocate
      for i=1:numel(st.bas)           % loop over basis set objects in unit cell
        b = st.bas{i}; ns = noff(i)+(1:b.Nf);
        % call generic strip discrep routine (which is overloaded for QP bases):
        Q(:,ns) = b.evalboundedstripdiscrep(st);
      end
    end % func

  end  % methods
  
  % ---------------------------------------------------------------------------
  methods(Static)
  end
end
