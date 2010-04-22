% QPSTRIP - create a quasi-periodic strip domain of (semi-)infinite y extent
%
% st = qpstrip(e, k) creates strip of (bi)infinite extent in the y direction
%   with translation/lattice vector e (should have non-zero x displacement).
%   k is the wavenumber in the strip. Left side wall passes through -e/2,
%   right side wall through e/2.
%
% st = qpstrip(e, k, opts) includes the following options:
%    opts.M = sets # quadr nodes on Sommerfeld contour for each FTy LP.
%    opts.seg,pm = uses a segment from -e/2 to e/2 in x-extent as a lower or
%        upper boundary for a semi-infinite strip. If goes left-to-right, it's
%        taken as a lower boundary, otherwise, upper.
%
%
% ISSUES:
%  * Will they exclude interior regions? Or allow upper/lower segment bdry?
%  * make better strip plot, etc
%
% See also DOMAIN
%
% (C) 2010 Alex Barnett
classdef qpstrip < handle & domain
  properties
    a                             % alpha QP phase factor
    e                             % horizontal lattice vector (complex #)
    r                             % reciprocal lattice vector (complex #)
    Lo, Ro                    % ftylayerpot origins for L, R walls 
    N                    % total # degrees of freedom in QP basis sets
    basnoff              % dof offsets for QP basis sets (numeric list)
    buffer               % 0,1: simple discrepancy or 3-unit-cell discrepancy
    up                   % +1 when seg is lower bdry, -1 if upper, 0 if no seg
    % Also seg, pm : L-to-R or R-to-L boundary segment & sense (optional)
  end
  
  methods
    function st = qpstrip(e, k, o) % ........................... constructor
      if nargin<3, o = []; end
      if real(e)<=0, error('lattice vector must have x component > 0!'); end
      st = st@domain();           % use R^2 domain
      st.perim = Inf; st.exterior = nan;   % since walls infinite extent in y
      st.Lo = -e/2; st.Ro = e/2;
      st.e = e;
      st.buffer = 0;
      st.a = 1;  % default (don't need st.setbloch method)
      st.k = k;
      st.r = 2*pi/st.e;  % works for real e at least
      if isfield(o,'seg'), st.seg = o.seg; st.pm = o.pm;
        tiny = 1e-14;        % semi-infinite in y extent...
        s1 = st.seg.eloc((o.pm+3)/2); s2 = st.seg.eloc(3-(o.pm+3)/2); % sta,end
        if abs(real(s1)-st.Lo)<tiny, st.up = -1;
          if abs(real(s2)-st.Ro)>tiny,
            error('seg end does not connect to R wall!'); end
        elseif abs(real(s1)-st.Ro)<tiny, st.up = +1;
          if abs(real(s2)-st.Lo)>tiny,
            error('set end does not connect to L wall!'); end
        else error('seg start connects to neither L nor R wall!');
        end
      end
    end

    function setbloch(st, a)
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
    % ADDQPFTYLAYERPOTS - add SLP+DLP FT-y layer potentials to L in qpstrip
    %
    % All calling opts are passed to ftylayerpot, apart from opts.omega which
    % is inherited from the qpstrip domain.
      if ~isfield(o,'omega'), o.omega = st.k; end      % inherit wavenumber
      % compute distance from origin to nearest singularity given Bloch params:
      if ~isfield(o,'nearsing')
        o.nearsing=min(abs(sqrt(st.k^2-(log(st.a)/1i+(-100:100)*2*pi/st.e).^2))); end
      % append the two types of layerpot to existing bases...
      st.bas = {st.bas{:}, ftylayerpot(st.Lo, 'd', o), ...
                ftylayerpot(st.Lo, 's', o)};
      for i=0:1, st.bas{end-i}.doms = st; end  % make strip the affected domain
      st.setupbasisdofs;
    end
    
    function requadrature(uc, N) % ....... overloads domain.requadrature
    % REQUADRATURE - reset # quadrature pts for Sommerfeld type FTy layerpots
    %
    % make this handle more general bases in the qpstrip...
      if nargin<2, N = []; end       % uses default
      if ~isa(uc.bas{1}, 'ftylayerpot')
        error('expecting bases 1,2 to be ftylayerpots...'); end
      uc.bas{1}.requadrature(N); uc.bas{2}.requadrature(N);
      uc.setupbasisdofs;
    end
    
    function i = inside(s, p) % ...................... overloaded from domain
    % INSIDE - compute if pointset locations fall inside (semi-bounded) strip
    %
    % Issues:
    % * currently hacked so only horiz segs handles correctly.
    % * need to fix up to handle excluded regions ? Think about it
    %
    % See also: domain.INSIDE
      i = real(p)>=s.Lo & real(p)<=s.Ro;
      if ~isempty(s.seg)
        yo = imag(s.seg.eloc(1));      % hack assuming all segs horizontal:
        i = i & ((imag(p)-yo)*s.up>0); % so, whether above or below the yo line
      end
      i = reshape(i, size(p));
    end
    
    function h = plot(s, varargin) % .................... plot it
    % PLOT - plot geometry of a qpstrip domain
      h = [plot@domain(s, varargin{:}); ...
           vline(s.Lo, 'k-', 'L'); vline(s.Ro, 'k-', 'R')]; %do as domain
    end
    
    function Q = evalbasesdiscrep(st, opts) % .................... Q
    % EVALBASESDISCREP - fill Q matrix mapping UC bases coeffs to discrepancy
    %
    %  This stacks matrices from basis.evalftystripdiscrep
    %
    %  Q = EVALBASESDISCREP(uc) returns matrix of basis function evaluations
    %   given basis set dofs as stored in unit cell.
    %
    % Based on: qpunitcell.evalbasesdiscrep
      if nargin<2, opts = []; end
      N = st.N; noff = st.basnoff;     % # basis dofs
      M = 2*st.bas{1}.Nf;              % # discrep dofs
      Q = zeros(M,N);                          % preallocate
      for i=1:numel(st.bas)           % loop over basis set objects in unit cell
        b = st.bas{i}; ns = noff(i)+(1:b.Nf);
        if ~isa(b,'rpwbasis')
          Q(:,ns) = b.evalftystripdiscrep(st); % currently only for ftylayerpots
        else
          Q(:,ns) = zeros(M, numel(ns)); % PW (delta in k) no effect on contour
        end
      end
    end % func

  end  % methods
  
  % ---------------------------------------------------------------------------
  methods(Static)
  end
end
