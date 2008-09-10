classdef layerpot < handle & basis

% LAYERPOT - create a layer potential basis set on a segment
%
%  b = LAYERPOT(seg, a, k, opts) where a = 'single' or 'double' creates a layer
%   potential basis object with density on segment seg, with wavenumber k
%   (which may be [] in which k is not defined), and options:
%    opts.real: if true, real valued (Y_0), otherwise complex (H_0 outgoing).
%    opts.fast: if 0 use matlab Hankel, 1 use fast fortran Hankel, 2 faster.
%   If a is instead a 1-by-2 array, then a mixture of a(1) single plus a(2)
%   double is created, useful for Brakhage, Werner, Leis & Panich type
%   representations.
%
%   Note that the discretization of the layerpot is given by that of the seg,
%    apart from periodic segments where new quadrature weights may be used.
%
%  To do: real=1 case? don't bother.

  properties
    real                            % true if fund sol is Y_0, false for H_0^1
    seg                             % handle of segment on which density sits
    a                               % 1-by-2, mixture weights of SLP and DLP
    quad                            % quadrature option (opts.quad in S,D,T)
    ord                             % quadrature order (opts.ord in S,D,T)
    fast                            % Hankel function evaluation method (0,1,2)
  end

  methods
    function b = layerpot(seg, a, k, opts) %........................ constructor
      if nargin<4, opts = []; end
      if nargin>2 & ~isempty(k), b.k = k; end
      if ~isfield(opts, 'fast'), opts.fast = 2; end
      % TODO: speed up the following check which happens every time created...
      if opts.fast==2 & exist('@utils/greengardrokhlinhank106.mexglx')~=3
        opts.fast = 1;               % downgrade the speed if 106 not available
      end
      b.fast = opts.fast;
      if ~isfield(opts, 'real'), opts.real = 0; end
      b.real = opts.real;
      if isfield(opts, 'quad'), b.quad = opts.quad; end  % quad=[] is default
      if isfield(opts, 'ord'), b.ord = opts.ord; end     % ord=[] is default
      if ~isa(seg, 'segment'), error('seg must be a segment object!'); end
      b.seg = seg;
      b.updateNf;                 % count the # dofs, is # quadr pts on seg
      b.N = b.Nf;
      if ~isnumeric(a)
        switch a
         case {'single', 'S', 's', 'SLP'}
          a = [1 0];
         case {'double', 'D', 'd', 'DLP'}
          a = [0 1];
         otherwise
          error('cannot interpret non-numeric input argument a!');
        end
      end
      if size(a)~=[1 2], error('a argument is not size 1-by-2!'); end
      b.a = a;
    end % func
    
    function updateNf(b)  % ............... Nf property reads segment # quadr
      b.Nf = numel(b.seg.x);
    end
    
    function [A Ax Ay] = eval(b, p, o) % .........basis evaluator at points p
% EVAL - evaluate layer potential on pointset or on a segment, with jump rel
%
%  A = EVAL(bas, p, opts) returns a matrix mapping degrees of freedom in the
%   discretization of the layer density to the values on p the target
%   pointset or segment.
%  (Includes source quadr wei, ie dofs are density values, as always for LP)
%
%  [A An] = EVAL(bas, p, opts) also returns normal derivatives using normals
%   in p.
%
%  [A Ax Ay] = EVAL(bas, p, opts) instead returns x- and y-derivatives,
%   ignoring the normals in p.
%
%   The optional argument opts may contain the following fields:
%    opts.dom: domain in which evaluation is performed. In the case of
%     p being the segment on which the LP density sits, opts.dom must be
%     defined, since it resolves which sign the jump relation has.
%    opts.fast: if 0 use matlab Hankel, 1 use fast fortran Hankel, 2 even faster
      if nargin<3, o = []; end
      if ~isfield(o, 'fast'), o.fast = b.fast; end    % default given in b obj
      b.updateNf;
      if ~isempty(b.quad), o.quad = b.quad; end % pass quadr type to S,D,T eval
      if ~isempty(b.ord), o.ord = b.ord; end % pass quadr order to S,D,T eval
      self = (p==b.seg); % if true, local eval on segment which carries density
      if self
        if isempty(o) | ~isfield(o, 'dom') | (o.dom~=p.dom{1} & o.dom~=p.dom{2})
          error('opts.dom must be a domain connected to segment p, since self eval with jump relation!');
        end
        % tell the jump relation from which side we're taking the limit...
        approachside = -1;           % - sign since normals point away from dom
        if o.dom==p.dom{1}, approachside = +1; end  % instead dom on + normal
        p = [];              % tell S, D, etc to use self-interaction
      end
      k = b.k;                          % could look up from b.doms someday...
      % precompute the displacement and distance matrices... used throughout
      N = numel(b.seg.x);                                % # src pts
      t = p; if self, t = b.seg; end                     % use src if self
      M = numel(t.x);                                    % # target pts
      o.displ = repmat(t.x, [1 N]) - repmat(b.seg.x.', [M 1]); % C-# displ mat
      o.rdist = abs(o.displ);
      if k>0 & o.fast % ...precompute Hankel H0,H1 vals in Sker and Dker_noang
        if o.fast==1
          [o.Sker o.Dker_noang] = utils.greengardrokhlinhank103(k*o.rdist);
        elseif o.fast==2
          [o.Sker o.Dker_noang] = utils.greengardrokhlinhank106(k*o.rdist);
        end
        o.Sker = (1i/4) * o.Sker;
        o.Dker_noang = (1i*k/4) * o.Dker_noang;
      end

      if nargout==1 %------------------------------- values only
        if b.a(2)==0             % only SLP
          A = layerpot.S(k, b.seg, p, o);
          if b.a(1)~=1, A = A * b.a(1); end 
        elseif b.a(1)==0         % only DLP
          A = layerpot.D(k, b.seg, p, o);
          if b.a(2)~=1, A = A * b.a(2); end 
        else                     % mixture of SLP + DLP
          A = b.a(1) * layerpot.S(k, b.seg, p, o) + ...
              b.a(2) * layerpot.D(k, b.seg, p, o);
        end
        if self & b.a(2)~=0   % NOTE should speed up by writing diag vals only:
          A = A + approachside * b.a(2) * eye(size(A)) / 2;  % DLP val jump
        end
        
      elseif nargout==2 %------------------------------- values + normal derivs
        if b.a(2)==0             % only SLP
          A = layerpot.S(k, b.seg, p, o);
          o.derivSLP = 1;
          Ax = layerpot.D(k, b.seg, p, o);
          if b.a(1)~=1, A = A * b.a(1); Ax = Ax * b.a(1); end 
        elseif b.a(1)==0         % only DLP
          A = layerpot.D(k, b.seg, p, o);
          Ax = layerpot.T(k, b.seg, p, o);
          if b.a(2)~=1, A = A * b.a(2); Ax = Ax * b.a(2); end 
        else                     % mixture of SLP + DLP
          [A o.Sker] = layerpot.S(k, b.seg, p, o);
          [AD o.Dker_noang] = layerpot.D(k, b.seg, p, o);
          A = b.a(1) * A + b.a(2) * AD;
          clear AD
          o.derivSLP = 1;
          [Ax o.Dker_noang] = layerpot.D(b.k, b.seg, p, o);
          Dn = layerpot.T(b.k, b.seg, p, o);
          Ax = b.a(1) * Ax + b.a(2) * Dn;
        end
        if self & b.a(2)~=0
          A = A + approachside * b.a(2) * eye(size(A)) / 2;   % DLP val jump
        end
        if self & b.a(1)~=0
          Ax = Ax - approachside * b.a(1) * eye(size(A)) / 2; % SLP deriv jump
        end
        
      else % ------------------------------------- values, x- and y-derivs
        if b.a(2)==0             % only SLP
          A = layerpot.S(k, b.seg, p, o);
          o.derivSLP = 1;
          nx_copy = p.nx; p.nx = ones(size(p.nx)); % overwrite pointset normals
          [Ax o.Dker_noang] = layerpot.D(k, b.seg, p, o);
          p.nx = 1i*ones(size(p.nx));
          Ay = layerpot.D(k, b.seg, p, o); % reuses H1 (Dker) mat
          p.nx = nx_copy;                          % restore pointset normals
          if b.a(1)~=1, A = A * b.a(1); Ax = Ax * b.a(1); Ay = Ay * b.a(1); end 
       elseif b.a(1)==0         % only DLP
          [A o.Dker_noang] = layerpot.D(k, b.seg, p, o);
          nx_copy = p.nx; p.nx = ones(size(p.nx)); % overwrite pointset normals
          [Ax o.Sker] = layerpot.T(k, b.seg, p, o);
          p.nx = 1i*ones(size(p.nx));
          Ay = layerpot.T(k, b.seg, p, o);
          p.nx = nx_copy;                          % restore pointset normals
          if b.a(2)~=1, A = A * b.a(2); Ax = Ax * b.a(2); Ay = Ay * b.a(2); end
        else                    % mixture of SLP + DLP
          [A o.Sker] = layerpot.S(k, b.seg, p, o);
          [AD o.Dker_noang] = layerpot.D(k, b.seg, p, o);
          A = b.a(1) * A + b.a(2) * AD;
          clear AD
          o.derivSLP = 1;
          nx_copy = p.nx; p.nx = ones(size(p.nx)); % overwrite pointset normals
          Ax = layerpot.D(k, b.seg, p, o);
          Ax = b.a(1) * Ax + b.a(2) * layerpot.T(k, b.seg, p, o);
          p.nx = 1i*ones(size(p.nx));
          Ay= layerpot.D(k, b.seg, p, o);
          Ay = b.a(1) * Ay + b.a(2) * layerpot.T(k, b.seg, p, o);
          p.nx = nx_copy;                          % restore pointset normals
        end
        if self
          error('LP self jump rels not yet implemented for Ax Ay grad vec!')
        end
        
      end
    end % func
    
    function showgeom(b, opts) % .................. crude show discr pts of seg
      b.seg.plot;
    end
    
  end % methods
  
  methods (Static) %...................... the beef: kernel evaluation routines
    [A Sker] = S(k, s, t, o)                      % Phi kernel matrix
    [A Dker_noang cosker] = D(k, s, t, o)         % dPhi/dny or dPhi/dny matrix
    [A Sker Dker_noang] = T(k, s, t, o)           % d^2Phi/dnxdny matrix
  end % methods
end