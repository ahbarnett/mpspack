% QPSTRIP - create a quasi-periodic strip domain of (semi-)infinite y extent
%
% st = qpstrip(e, k) creates strip of (bi)infinite extent in the y direction
%   with translation/lattice vector e (should have non-zero x displacement).
%   k is the wavenumber in the strip. Left side wall passes through origin,
%   right side wall through e.
%
% ISSUES:
%  * Will they exclude interior regions? Or allow upper/lower segment bdry?
%  * make better strip plot, etc
%
% See also DOMAIN
%
% (C) Alex Barnett
classdef qpstrip < handle & domain
  properties
    a                             % alpha QP phase factor
    e                             % horizontal lattice vector (complex #)
    r                             % reciprocal lattice vector (complex #)
    Lo, Ro                    % ftylayerpot origins for L, R walls 
    N                    % total # degrees of freedom in QP basis sets
    basnoff              % dof offsets for QP basis sets (numeric list)
    buffer               % 0,1: simple discrepancy or 3-unit-cell discrepancy
  end
  
  methods
    function st = qpstrip(e, k, opts) % ........................... constructor
      if real(e)<=0, error('lattice vector must have x component > 0!'); end
      st = st@domain();           % use R^2 domain
      st.perim = Inf; st.exterior = nan;   % since walls infinite extent in y
      st.Lo = 0; st.Ro = e;
      st.e = e;
      st.buffer = 0;
      st.a = 1;  % default (don't need st.setbloch method)
      st.k = k;
    end

    function setbloch(st, a)
      st.a = a;          % TODO: need to add recip
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

    function addqpftylayerpots(st, varargin) % ...............................
    % ADDQPFTYLAYERPOTS - add SLP+DLP FT-y layer potentials to L in qpstrip
      if st.buffer==0
        st.bas = {st.bas{:}, ftylayerpot(st.Lo, 'd', varargin{:}), ...
                  ftylayerpot(st.Lo, 's', varargin{:})};
      else
        error('buffer>0 not implemented!');
      end
      for i=0:1, st.bas{end-i}.doms = st; end  % make strip the affected domain
      st.setupbasisdofs;
    end
    
   function requadrature(uc, N) % ....... overloads domain.requadrature
    % REQUADRATURE - reset # quadrature pts for Sommerfeld type FTy layerpots
      if nargin<2, N = []; end       % uses default
      uc.bas{1}.requadrature(N); uc.bas{2}.requadrature(N);
      uc.setupbasisdofs;
    end
    
    function i = inside(s, p) % ...................... overloaded from domain
    % need to fix up to handle excluded regions ? Think about it
      i = real(p)>=0 & real(p)<=real(s.e);
      i = reshape(i, size(p));
    end
    
    function h = plot(s, varargin) % .................... plot it
    % PLOT - plot geometry of a qpstrip domain
      h = [plot@domain(s, varargin{:}); ...
           vline(s.Lo, 'r-', 'L'); vline(s.Ro, 'r-', 'R')]; %do as domain
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
      % inititialize struct cell array to be used to store unphased Q data...
      for i=1:numel(st.bas)           % loop over basis set objects in unit cell
        b = st.bas{i}; ns = noff(i)+(1:b.Nf);
        % call generic discrep routine (which is overloaded for QP bases):
        Q(:,ns) = b.evalftystripdiscrep(st);
      end
    end % func

  end  % methods
  
  % ---------------------------------------------------------------------------
  methods(Static)
  end
end
