classdef layerpot < handle & basis

% LAYERPOT - create a layer potential basis set on a segment
%
%  b = LAYERPOT(seg, a, k, opts) where a = 'single' or 'double' creates a layer
%   potential basis object with density on segment seg, with wavenumber k
%   (which may be [] in which k is not defined), and options:
%   opts.real: if true, real valued (Y_0), otherwise complex (H_0 outgoing).
%   If a is instead a 1-by-2 array, then a mixture of a(1) single plus a(2)
%   double is created, useful for Brakhage, Werner, Leis & Panich type
%   representations.
%
%   Note that the discretization of the layerpot is given by that of the seg,
%    apart from periodic segments where new quadrature weights may be used.
%
%  To do: real=1 case?

  properties
    real                            % true if fund sol is Y_0, false for H_0^1
    seg                             % handle of segment on which density sits
    a                               % 1-by-2, mixture weights of SLP and DLP
  end

  methods
    function b = layerpot(seg, a, k, opts) %........................ constructor
      if nargin<4, opts = []; end
      if nargin>2 & ~isempty(k), b.k = k; end
      b.N = numel(seg.x);
      b.Nf = b.N;
      if ~isfield(opts, 'real'), opts.real = 0; end
      b.real = opts.real;
      if ~isa(seg, 'segment'), error('seg must be a segment object!'); end
      b.seg = seg;
      switch a
        case {'single', 'S', 's', 'SLP'}
         a = [1 0];
        case {'double', 'D', 'd', 'DLP'}
         a = [0 1];
      end
      b.a = a;
    end % func
    
    function [A Ax Ay] = eval(b, p, o) % .........basis evaluator at points p
    % EVAL - evaluate layer potential on pointset or on a segment, with jump rel
    %
    %  A = EVAL(bas, p, opts)           p is the target pointset or segment
    %  [A An] = EVAL(bas, p, opts)
    %  [A Ax Ay] = EVAL(bas, p, opts)
    %
    %   The optional argument opts may contain the following fields:
    %    opts.layerpotside = +1 or -1: determines from which side of a segment
    %     the limit is approached, for jump relation. +1 is the normal side.
      self = (p==b.seg); % if true, local eval on segment carrying density
      if self & o.layerpotside~=1 & o.layerpotside~=-1
        error('opts.layerpotside must be +1 or -1 for self-interaction');
      end
      if self, p = []; end              % tell S, D, etc to use self-interaction
      if nargout==1 %------------------------------- values only
        if b.a(2)==0             % only SLP
          A = b.a(1) * layerpot.S(b.k, b.seg, p, o);
        elseif b.a(1)==0         % only DLP
          A = b.a(2) * layerpot.D(b.k, b.seg, p, o);
        else                     % mixture of SLP + DLP
          A = b.a(1) * layerpot.S(b.k, b.seg, p, o) + ...
              b.a(2) * layerpot.D(b.k, b.seg, p, o);
        end
        if self & b.a(2)~=0
          A = A + o.layerpotside * b.a(2) * eye(size(A)) / 2;  % DLP val jump
        end
        
      elseif nargout==2 %------------------------------- values + normal derivs
        
        
        
      else
        
        
      end
    end % func
  end % methods
  
  methods (Static)  
    function A = S(k, s, t, o)
    % S - double layer potential discretization matrix for density on a segment
    %
    %  S = S(k, s) where s is a segment returns the matrix discretization of
    %   the SLP integral operator with density on the segment,
    %      u(x) = int_s Phi(x,y) sigma(y) ds(y)
    %   where Phi is the fundamental solution
    %      Phi(x,y) = (i/4) H_0^{(1)}(k|x-y|)    for k>0
    %                 (1/2pi) log(1/|x-y|)       for k=0
    %   for x each of the points on the self-same segment s. No jump relations
    %   are included. If the segment is closed then periodic spectral
    %   quadrature may be used (see options below). 
    %
    %  S = S(k, s, t) where t is a pointset (any object with fields t.x)
    %   uses s as the source segment as above, but target points in t.
    %   It is assumed t is not s.
    %
    %  S = S(k, s, [], opts) or S = S(k, s, t, opts) does the above two choices
    %   but with options struct opts including the following:
    %    opts.quad = 'k' (Kapur-Rokhlin), 'm' (Martensen-Kussmaul spectral)
    %                periodic quadrature rules, used only if s.qtype is 'p'
    %    opts.ord = 2,6,10. controls order of Kapur-Rokhlin rule.
    %    opts.Aval = quad-unweighted value matrix of fund-sols (prevents
    %                recomputation of r or fundsol values). Non-self case only!
    %
    %  Adapted from leslie/pc2d/slp_matrix.m, barnett 7/31/08
    %  Tested by routine: testlpquad.m
      if isempty(k) | isnan(k), error('SLP: k must be a number'); end
      if nargin<4, o = []; end
      if nargin<3, t = []; end
      if ~isfield(o, 'quad'), o.quad='m'; end; % default self periodic quadr
      self = isempty(t);               % self-interact: potentially sing kernel
      N = numel(s.x);                  % # src pts
      if self, t = s; end              % use source as target pts (handle copy)
      M = numel(t.x);                  % # target pts
      sp = s.speed/2/pi; % note: 2pi converts to speed wrt s in [0,2pi] @ src pt
      needAval = self | ~isfield(o, 'Aval'); % true if must compute kernel vals
      if needAval
        d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]); % C-# displacements mat
        r = abs(d);                                    % dist matrix R^{MxN}
      end
      if self, r(diagind(r)) = 999; end   % dummy nonzero diag values

      if self % ........... source curve = target curve; can be singular kernel
  
        if s.qtype=='p' & o.quad=='k'  % Kapur-Rokhlin (kills diagonal values)
          A = utils.fundsol(r, k);              % Phi
          [s w] = quadr.kapurtrap(N+1, o.ord);  % Zydrunas-supplied K-R weights
          w = 2*pi * w;                 % change interval from [0,1) to [0,2pi)
          A = circulant(w(1:end-1)) .* A .* repmat(sp.', [M 1]); % speed
    
        elseif s.qtype=='p' & o.quad=='m' % Martensen-Kussmaul (Kress MCM 1991)
          if isempty(s.kappa)
            error('cant do Martensen-Kussmaul quadr without s.kappa available!')
          end
          A = utils.fundsol(r, k);                % Phi
          if k==0                         % Laplace SLP has ln sing
            S1 = -1/4/pi;                 % const M_1/2 of Kress w/out speed fac
            A = A - S1.*circulant(log(4*sin(pi*(0:N-1)/N).^2)); % A=D2=M_2/2 "
            A(diagind(A)) = -log(sp)/2/pi;        % diag vals propto curvature
          else
            S1 = triu(besselj(0,k*triu(r,1)),1);  % use symmetry (arg=0 is fast)
            S1 = -(1/4/pi)*(S1.'+S1);     % next fix it as if diag(r) were 0
            S1(diagind(S1)) = -(1/4/pi);  % S1=M_1/2 of Kress w/out speed fac
            A = A - S1.*circulant(log(4*sin(pi*(0:N-1)/N).^2)); % A=D2=M_2/2 "
            eulergamma = -psi(1);         % now set diag vals Kress M_2(t,t)/2
            A(diagind(A)) = 1i/4 - eulergamma/2/pi - log((k*sp).^2/4)/4/pi;
          end
          %if N==450, figure; imagesc(real(A)); colorbar; end % diag matches?
          A = (circulant(quadr.kress_Rjn(N/2)).*S1 + A*(2*pi/N)) .* ...
              repmat(sp.', [M 1]);  

        else       % ------ self-interacts, but no special quadr, just use seg's
          % Use the crude approximation of kappa for diag, usual s.w off-diag...
          A = utils.fundsol(r, k);
          A = A .* repmat(s.w, [M 1]);  % use segment usual quadrature weights
          fprintf('warning: SLP crude self-quadr will be awful, no diag!\n')
          A(diagind(A)) = 0;
        end
      
      else % ............................ distant target curve, so smooth kernel

        if needAval
          A = utils.fundsol(r, k);                      
          Aval = A;
        else
          A = o.Aval;
        end
        A = A .* repmat(s.w, [M 1]);       % use segment quadrature weights
      end
    end % func

    function A = D(k, s, t, o)
    % D - double layer potential discretization matrix for density on a segment
    %
    %  D = D(k, s) where s is a segment returns the matrix discretization of
    %   the DLP integral operator with density on the segment,
    %      u(x) = int_s partial_n(y) Phi(x,y) tau(y) ds(y)
    %   where Phi is the fundamental solution
    %      Phi(x,y) = (i/4) H_0^{(1)}(k|x-y|)    for k>0
    %                 (1/2pi) log(1/|x-y|)       for k=0
    %   for x each of the points on the self-same segment s. No jump relations
    %   are included. If the segment is closed then periodic spectral
    %   quadrature may be used (see options below). 
    %
    %  D = D(k, s, t) where t is a pointset (any object with fields t.x)
    %   uses s as the source segment as above, but target points in t.
    %   It is assumed t is not s.
    %
    %  D = D(k, s, [], opts) or D = D(k, s, t, opts) does the above two choices
    %   but with options struct opts including the following:
    %    opts.quad = 'k' (Kapur-Rokhlin), 'm' (Martensen-Kussmaul spectral)
    %                periodic quadrature rules, used only if s.qtype is 'p'
    %    opts.ord = 2,6,10. controls order of Kapur-Rokhlin rule.
    %
    %  Adapted from leslie/pc2d/dlp_matrix.m, barnett 7/31/08
    %  Tested by routine: testlpquad.m
    %
    %  To Do: * put option in to return D^T instead, switching s.nx to t.nx
    %         * make it bypass hankel computation if pass in opts.Hval matrix!
      if isempty(k) | isnan(k), error('DLP: k must be a number'); end
      if nargin<4, o = []; end
      if nargin<3, t = []; end
      if ~isfield(o, 'quad'), o.quad='m'; end; % default self periodic quadr
      self = isempty(t);               % self-interact: potentially sing kernel
      N = numel(s.x);                  % # src pts
      if self, t = s; end              % use source as target pts (handle copy)
      M = numel(t.x);                  % # target pts
      sp = s.speed/2/pi; % note: 2pi converts to speed wrt s in [0,2pi] @ src pt
      d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]); % C-# displacements matrix
      nx = repmat(s.nx.', [M 1]);         % identical rows given by src normals
      r = abs(d);                         % dist matrix R^{MxN}
      if self, r(diagind(r)) = 999; end   % dummy nonzero diag values
      cosphi = real(conj(nx).*d) ./ r;    % dot prod <normal, displacement>

      if self % ........... source curve = target curve; can be singular kernel
  
        if s.qtype=='p' & o.quad=='k' % Kapur-Rokhlin (kills diagonal values)
          A = utils.fundsol_deriv(r, cosphi, k);      % (src)normal-deriv of Phi
          [s w] = quadr.kapurtrap(N+1, o.ord);  % Zydrunas-supplied K-R weights
          w = 2*pi * w;                 % change interval from [0,1) to [0,2pi)
          A = circulant(w(1:end-1)) .* A .* repmat(sp.', [M 1]); % speed
    
        elseif s.qtype=='p' & o.quad=='m' % Martensen-Kussmaul (Kress MCM 1991)
          if isempty(s.kappa)
            error('cant do Martensen-Kussmaul quadr without s.kappa available!')
          end
          if k==0                   % for k=0 DLP analytic on analytic curve
            A = utils.fundsol_deriv(r, cosphi, k);  % (src)normal-deriv of Phi
            A(diagind(A)) = -s.kappa/(4*pi);        % diag vals propto curvature
            A = (2*pi/N) * A .* repmat(sp.', [M 1]);   % speed quad weights
          else                      % for k>0 has x^2 ln x sing, Kress handles
            A = utils.fundsol_deriv(r, cosphi, k);     % diag is garbage
            D1 = triu(besselj(1,k*triu(r,1)),1); % use symmetry (arg=0 is fast)
            D1 = -(k/4/pi)*cosphi.*(D1.'+D1);  % L_1/2 of Kress w/out speed fac
            A = A - D1.*circulant(log(4*sin(pi*(0:N-1)/N).^2)); % A=D2=L_2/2 "
            A(diagind(A)) = -s.kappa/(4*pi);   % L_2(t,t)/2, same as for k=0
            % speed factors: diag matrix mult from right...
            A = (circulant(quadr.kress_Rjn(N/2)).*D1 + (2*pi/N)*A) .* ...
                repmat(sp.', [M 1]);
          end

        else       % ------ self-interacts, but no special quadr, just use seg's
          % Use the crude approximation of kappa for diag, usual s.w off-diag...
          A = utils.fundsol_deriv(r, cosphi, k); % src-normal-deriv of fund soln
          A = A .* repmat(s.w, [M 1]);  % use segment usual quadrature weights
          if isempty(s.kappa)
            fprintf('warning: DLP crude self-quadr has no s.kappa so will be awful!\n')
            A(diagind(A)) = 0;
          else
            fprintf('warning: DLP crude self-quadr, using s.kappa for diag\n')
            A(diagind(A)) = -s.kappa/(4*pi);       % diag vals propto curvature
          end
        end
      
      else % ............................ distant target curve, so smooth kernel

        A = utils.fundsol_deriv(r, cosphi, k);   % src-normal-deriv of fund soln
        A = A .* repmat(s.w, [M 1]);       % use segment quadrature weights
      end
    end % func
  end % methods
end