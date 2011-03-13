classdef layerpot < handle & basis

% LAYERPOT - create a layer potential basis set on a segment
%
%  b = LAYERPOT(seg, a, opts) where a = 'single' or 'double' creates a layer
%   potential basis object with density on segment seg, with options:
%    opts.real: if true, real valued (Y_0), otherwise complex (H_0 outgoing).
%    opts.fast: if 0 use matlab Hankel, 1 use fast fortran Hankel, 2 faster.
%   If a is instead a 1-by-2 array, then a mixture of a(1) single plus a(2)
%   double is created, useful for Brakhage, Werner, Leis & Panich type
%   combined representations. As usual for basis sets, the wavenumber k is
%   determined by that of the affected domain.
%
%   Note that the discretization of the layerpot is given by that of the seg,
%    apart from periodic segments where new quadrature weights may be used.
%
%
% See also: DOMAIN/ADDLAYERPOT

% Copyright (C) 2008-2011, Alex Barnett and Timo Betcke
  properties
    real                            % true if fund sol is Y_0, false for H_0^1
    seg                             % handle of segment on which density sits
    a                               % 1-by-2, mixture weights of SLP and DLP
    quad                            % quadrature option (opts.quad in S,D,T)
    ord                             % quadrature order (opts.ord in S,D,T)
    fast                            % Hankel function evaluation method (0,1,2)
    Jfilter                         % (optional) Jfilter structure
    self                            % self-int matrix struct, w/ S, D, DT, T
  end

  methods
    function b = layerpot(seg, a, opts) %........................ constructor
      if nargin<3, opts = []; end
      if ~isfield(opts, 'fast'), opts.fast = 2; end
      % TODO: speed up the following check which happens every time created...
      if opts.fast==2 & exist('greengardrokhlinhank106')~=3
        opts.fast = 1;               % downgrade the speed if 106 not available
      end
      if opts.fast==1 & exist('greengardrokhlinhank103')~=3
        opts.fast = 0;               % downgrade the speed if 103 not available
      end
      b.fast = opts.fast;
      b.HFMMable = (exist('hfmm2d')==3);  % is Helmholtz FMM MEX available?
      if ~isfield(opts, 'real'), opts.real = 0; end
      b.real = opts.real;
      if isfield(opts, 'quad'), b.quad = opts.quad; end  % quad=[] is default
      if isfield(opts, 'ord'), b.ord = opts.ord; end     % ord=[] is default
      if ~isfield(opts,'nmultiplier'), opts.nmultiplier=1; end
      b.nmultiplier = opts.nmultiplier;
      if ~isa(seg, 'segment'), error('seg must be a segment object!'); end
      b.seg = seg;
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
    
    function Nf = Nf(b)  % ............... Nf method reads segment # quadr
      Nf = numel(b.seg.x);
      % was this too slow, some reason why Nf was a property ???
    end
    
    function updateN(b,N) % ................ overloaded in some basis objects
    % UPDATEN - Change N and requadrature segment in proportion to an overall N
      b.N = ceil(N * b.nmultiplier / 2) * 2;  % ensures even #, for Kress quad
      b.seg.requadrature(b.N);
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
%    opts.fast: if 0 use matlab Hankel, 1,2 use Rokhlin's fortran Hankels.
%    opts.Jfilter: (advanced feature) evaluate via J-expansion instead.
%     Expects structure with:
%                  Jfilter.M = max order of local expansion
%       (optional) Jfilter.rescale_rad = radius to rescale J-exp to be O(1) at
%       (optional) Jfilter.origin = center of J-expansion
%
%  To do: make it handle a segment list, with correct accounting of self-eval.
%
% See also: LAYERPOT.FMMEVAL
%
% (c) Alex Barnett 2008-2011
      if nargin<3, o = []; end
      if ~isfield(o, 'fast'), o.fast = b.fast; end    % default given in b obj
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
        if ~isempty(b.self), o.self = b.self; end    % pass in self-int mats
      end
      k = b.k(o);               % method gets k from affected domain (o.dom)
      %fprintf('o.self=%d, k=%g, type [%g %g]\n',self,k,b.a(1),b.a(2))
      
      % opts takes precedence, but otherwise use Jfilter in basis properties:
      if ~isempty(b.Jfilter) & ~isfield(o, 'Jfilter'), o.Jfilter=b.Jfilter; end
      if isfield(o, 'Jfilter') & ~isempty(p)  % ==== J-expansion filter, nonself
        
        if isempty(p), p=b.seg;      % experiment to do self-eval w/ J-filter
          warning 'layerpot.eval J-filter cannot do self-eval!'% probably true!
        end
        if isfield(o.Jfilter, 'rescale_rad')
          Jopts.rescale_rad = o.Jfilter.rescale_rad; % same resc for S2L as J
        end
        xo = 0; if isfield(o.Jfilter, 'origin'), xo=o.Jfilter.origin; end
        Jopts.real = 0; Jopts.besselcode = 'm';  % GSL needed for small J vals
        Jexp = regfbbasis(xo, o.Jfilter.M, Jopts); % make new local exp basis
        Jexp.doms = b.doms;      % FIX to include o.dom ???
               
        % fill S2L matrix (layer potential quadr pts to top level J-expansion):
        if b.a(2)==0             % only SLP
          S2L = layerpot.localfromSLP(b.seg, Jexp);  % k, rescale_rad from Jexp
          if b.a(1)~=1, S2L = S2L * b.a(1); end 
        elseif b.a(1)==0         % only DLP
          S2L = layerpot.localfromDLP(b.seg, Jexp);  % k, rescale_rad from Jexp
          if b.a(2)~=1, S2L = S2L * b.a(2); end
           S2L = layerpot.localfromDLP(b.seg, Jexp);
        else                     % mixture of SLP + DLP
          S2L = b.a(1) * layerpot.localfromSLP(b.seg, Jexp) + ...
                b.a(2) * layerpot.localfromDLP(b.seg, Jexp);
        end

        %max(abs(S2L(:)))
        %figure; imagesc(abs(S2L)); colorbar
        
        if nargout==1  % evaluate the local J-expansion, with correct # out args
          AJ = Jexp.eval(p);
          A = AJ * S2L;    % note O(N^3) .... don't know how to avoid yet ?
        elseif nargout==2
          [AJ AJn] = Jexp.eval(p);
          A = AJ * S2L; Ax = AJn * S2L;
        else
          [AJ AJx AJy] = Jexp.eval(p);
          A = AJ * S2L; Ax = AJx * S2L; Ay = AJy * S2L;
        end
      
      else  % ================= use usual layer-potential evaluations
      
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
        if self & b.a(2)~=0   % Jump Relation
          A(diagind(A)) = A(diagind(A)) + approachside*b.a(2)/2; % DLP val jump
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
        if self & b.a(2)~=0    % Jump Relations: DLP value jump
          A(diagind(A)) = A(diagind(A)) + approachside*b.a(2)/2;
        end
        if self & b.a(1)~=0    % SLP normal-derivative jump
          Ax(diagind(Ax)) = Ax(diagind(Ax)) - approachside*b.a(1)/2;
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
        
      end % ------ end nargout switch
      end % ================= end Jfilter switch
    end % func

    function [u u1 u2] = evalFMM(b, co, p, o)
      % EVALFMM - FMM evaluate layer potential (SLP+DLP) basis at set of points
      %
      % u = evalFMM(b, p, c) where b is a layerpot object, p is a pointset
      %  object containing M points, and c is a coefficient vector, returns
      %  u, the layer potential values at each target point m=1...M as a column
      %  vector. (Note p may also be a column vector of points in the C plane.)
      %  The Helmholtz fast multipole method of Gimbutas-Greengard is called.
      %
      % [u un] = evalFMM(b, p, c) also returns gradient of potential in the
      %  normal directions contained in the pointset p.nx.
      %
      % [u ux uy] = evalFMM(b, p, c) also returns gradient of potential in the
      %  x and y directions, ie its partials. (pointset normals not used).
      %
      % [u ...] = evalFMM(b, p, c, opts) controls options.
      %   The optional argument opts may contain the following fields:
      %    opts.dom: domain in which evaluation is performed. In the case of
      %     p being the segment on which the LP density sits, opts.dom must be
      %     defined, since it resolves which sign the jump relation has.
      %
      %  To do: make it handle a segment list target, with correct
      %         accounting of self-eval. Pass all non-self targs in a single
      %         FMM pass so can handle large numbers of segments efficiently
      %
      % See also: EVAL
      if nargin<4, o = []; end
      if isnumeric(p), x = p; else x = p.x; end    % allow point list target
      self = (p==b.seg); % if true, local eval on segment which carries density
      if self
        if isempty(o) | ~isfield(o, 'dom') | (o.dom~=p.dom{1} & o.dom~=p.dom{2})
          error('opts.dom must be a domain connected to segment p, since self eval with jump relation!');
        end
        % tell the jump relation from which side we're taking the limit...
        approachside = -1;           % - sign since normals point away from dom
        if o.dom==p.dom{1}, approachside = +1; end  % instead dom on + normal
      end
      k = b.k(o);               % method gets k from affected domain (o.dom)
      ifcharge = (b.a(1)~=0); ifdipole = (b.a(2)~=0);
      s = b.seg; N = numel(s.x); M = numel(x);  % NB s.x & s.nx are col vecs:
      source = [real(s.x) imag(s.x)].'; dipvec = [real(s.nx) imag(s.nx)].';
      target = [real(x) imag(x)].';    % since x is col vecs, but want 2-by-N
      iffldtarg = (nargout>1);
      iprec=4;                    % 12 digit precision - should be an opts
      sc = -1;                    % factor to convert to my i.H_0/4 scaling.
      charge = s.w .* co.';       % charge strengths = quadr weights * density

      if ~self                % non-self target (use segment's own quadrature)
        U = utils.hfmm2dparttarg(iprec,k,N,source,ifcharge,b.a(1)*charge,...
                  ifdipole,-b.a(2)*charge,dipvec,0,0,0,M,target,1,iffldtarg,0);
              % note the dipole strengths vector above has sign change!
        u = sc * U.pottarg.';  % as for mfsbasis, scale & convert to column vec
        if nargout==2          % field is supposedly -grad(potential) ...sign?
          u1 = sc * (real(p.nx).'.*U.fldtarg(1,:) + ...
                   imag(p.nx).'.*U.fldtarg(2,:)).';
        elseif nargout==3
          u1 = sc * U.fldtarg(1,:).'; u2 = sc * U.fldtarg(2,:).';
        end
      
      else                    % self-interaction
        if nargout==1
          u = utils.hfmm2dpart(iprec,k,N,source,ifcharge,b.a(1)*charge,...
                               ifdipole,-b.a(2)*charge,dipvec);
          u = sc * u.';   % scale & convert to col vec
          if self && b.a(2)~=0   % Jump Relation
            u = u + approachside*b.a(2)/2; % DLP val jump
          end
       % add in effect of quadratures...
       
       
       
       elseif nargout==2
          [u fld] = utils.hfmm2dpart(iprec,k,N,source,ifcharge,b.a(1)*charge,...
                                     ifdipole,-b.a(2)*charge,dipvec);
          u = sc * u.';   % scale & convert to col vec
          u1 = sc * (real(p.nx).'.*fld(1,:) + imag(p.nx).'.*fld(2,:)).';
          if self && b.a(2)~=0   % Jump Relations
            u = u + approachside*b.a(2)/2; % DLP val jump
          end
          if self && b.a(1)~=0
            u1 = u1 - approachside*b.a(1)/2; % SLP deriv jump
          end
        % add in effect of quadratures...
        
        
        
      else
          error('self interaction not implemented for returning [u ux uy]!');
        end
      end
    end
      
      
    function showgeom(b, opts) % .................. crude show discr pts of seg
      b.seg.plot;
    end
    
    % methods in separate files...
    [A Ax] = evalfty(b, f, o)                     % for 1.5D scatt, qpstrip
   
  end % methods
  
  methods (Static) %...................... the beef: kernel evaluation routines
    [A Sker] = S(k, s, t, o)                      % Phi kernel matrix
    [A Dker_noang cosker] = D(k, s, t, o)         % dPhi/dny or dPhi/dny matrix
    [A Sker Dker_noang] = T(k, s, t, o)           % d^2Phi/dnxdny matrix    
    A = localfromSLP(k, s, M, o)                  % Jfilter: S2L for SLP
    A = localfromDLP(k, s, M, o)                  % Jfilter: S2L for DLP
    A = fundsol(r, k)                             % to replace w/ utils.fundsol
    [B radderivs] = fundsol_deriv(r, cosphi, k, radderivs)
  end % methods
end
