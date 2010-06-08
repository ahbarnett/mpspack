% QPFTYLAYERPOT - set up L and R phased pair of FTy layer-pot bases on QP strip
%
% bas = QPFTYLAYERPOT(st, a) creates a pair of FTy layer potentials on L and R
%  walls of a qpstrip object st. a controls the amounts of SLP and DLP, as usual
%
% This is a simple class which basically duplicates the L wall layerpot onto
%  the R wall with Bloch phase from the underlying qpstrip object.
%
% See also: FTYLAYERPOT, QPSTRIP

% (C) 2010 Alex Barnett

classdef qpftylayerpot < handle & basis
  properties
    t                                            % qpstrip on which locs based
    lp                                           % layer-pot object on L seg
  end
  
  methods
    function b = qpftylayerpot(t, a, opts) % ................. constructor
      if nargin<3, opts = []; end
      b.t = t;
      b.lp = ftylayerpot(t.Lo, a, opts);
      b.lp.doms = t;
    end
    
    function Nf = Nf(b)  % ............... Nf method reads # quadr pts of L
      Nf = numel(b.lp.kj);
    end

    function requadrature(b, N, opts) % ................. requad
    % REQUADRATURE - make N-pt Sommerfeld contour quadrature for qpftlayerpot.
      if nargin<3, opts = []; end
      if nargin<2, error('must supply N for requadrature'); end
      requadrature(b.lp, N, opts);
    end

    function updateN(b,N) % ................ overloads from basis
    % UPDATEN - Change basis set # quadr pts in proportion to an overall N.
      updateN(b.lp,N);
    end

    function [A Ax Ay] = eval(b, p, o) % .........basis evaluator at points p
% EVAL - evaluate layer potential on pointset or on a segment, with jump rel
%
%  A = EVAL(bas, p, opts) returns a matrix mapping degrees of freedom in the
%   discretization of the Sommerfeld integral to the values on p the target
%   pointset or segment.
%
%  [A An] = EVAL(bas, p, opts) also returns normal derivatives using normals
%   in p.
%
%  [A Ax Ay] = EVAL(bas, p, opts) instead returns x- and y-derivatives,
%   ignoring the normals in p.
%
%  Bloch alpha is taken from the basis' underlying qpstrip object.
%
% See other options in: FTYLAYERPOT.EVAL
      if nargin<3, o = []; end
      pt = pointset(p.x - b.t.e, p.nx); % translated targets as if from R wall
      a = b.t.a;                        % current Bloch phase
      % use stuff from B filling in qpscatt.fillbcmatrix
      if nargout==1
        A = b.lp.eval(p, o) + a * b.lp.eval(pt, o);
      elseif nargout==2
        [A Ax] = b.lp.eval(p, o); [At Axt] = b.lp.eval(pt, o);
        A = A + a*At; Ax = Ax + a*Axt;  % note x really means n for this case
      else
        [A Ax Ay] = b.lp.eval(p, o); [At Axt Ayt] = b.lp.eval(pt, o);
        A = A + a*At; Ax = Ax + a*Axt; Ay = Ay + a*Ayt;
      end
    end
    
    function Q = evalftystripdiscrep(b, st) % ........... overloads from basis
    % EVALFTYSTRIPDISCREP - Q submatrix: QP FTy LP effect on FTy discrepancy 
    %
    % Q = EVALFTYSTRIPDISCREP(bas, st). bas is the basis set, st the strip
    %   object. Currently supports standard qpftylp's on L and copies on L+e,
    %   with buf=0,1... Q is a 2NxN matrix, with NxN diagonal subblocks, since
    %   transl ops are diag in FTy rep. Note, Bloch alpha is taken from st,
    %   rather than bas.t (the basis' underlying qpstrip object).
    % 
    % ISSUES: check the non-x-parallel e case (not urgent).
    %
    % * Obviously should make a multiply-by-Q routine since it's block diagonal!
      if ~isa(st, 'qpstrip'), error('st must be a qpstrip object!'); end
      N = b.Nf;
      om = b.k;                               % get omega from basis (ie domain)
      a = st.a;                               % strip's Bloch phase (neq pr.a?)
      d = real(st.e);                         % x-displacement across strip
      som = sqrt(om^2 - b.lp.kj.^2);          % Sommerfeld, row vec
      exf = exp(1i*som*d);                    % decay factors, row vec
      ba = b.lp.a;                            % S/DLP coeffs of FTy LPs
      Q = [diag((ba(1)*1i*(a-1/a)./som + ba(2)*(-a-1/a)).*exf/2); ...
           diag((ba(1)*(a+1/a) + ba(2)*1i*(a-1/a)*som).*exf/2)];
      i = diagind(Q); Q(i) = Q(i) + ba(2);   % jump relations add identities
      i = diagind(Q,N); Q(i) = Q(i) - ba(1);
    end
    
    function showgeom(b, opts) % .................. show location of L, R rays
      vline(real(b.lp.orig + [0 b.t.e]), 'k-');
    end
    
  end % methods
  
  methods (Static) %......................
  end % methods
end
