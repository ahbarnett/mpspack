function [A Sker] = S(k, s, t, o)
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
%               'a' (Alpert local correction)
%                periodic quadrature rules, used only if s.qtype is 'p';
%                any other does low-order non-periodic quad using segment's own.
%    opts.ord controls order of Kapur-Rokhlin (2,6,10) or Alpert rule (2,...,16)
%    opts.Sker = quad-unweighted kernel matrix of fund-sols (prevents
%                recomputation of r or fundsol values).
%    opts.rdist = matrix of source-target distances, prevents recalculation.
%    opts.close = distance below which adaptive quadrature is used to evaluate
%                distant targets (slow). If not present, never adaptive.
%    opts.closeacc = set relative tolerance for close eval (default 1e-12)
%    opts.self = self-interaction matrices, only self.S will be used.
%
% [S Sker] = S(...) also returns quad-unweighted kernel values matrix Sker.

% Copyright (C) 2008 - 2011, Alex Barnett and Timo Betcke

if isempty(k) | isnan(k), error('SLP: k must be a number'); return; end
if nargin<4, o = []; end
if nargin<3, t = []; end
if ~isfield(o, 'quad'), o.quad='m'; end; % default self periodic quadr
if ~isfield(o, 'ord'), o.ord=10; end;    % default quadr order
self = isempty(t);               % self-interact: potentially sing kernel
N = numel(s.x);                  % # src pts
if self, t = s; end              % use source as target pts (handle copy)
M = numel(t.x);                  % # target pts
sp = s.speed/2/pi; % note: 2pi converts to speed wrt s in [0,2pi] @ src pt
needA = ~isfield(o, 'Sker');     % true if must compute kernel vals
if isfield(o, 'rdist'), r = o.rdist; end
if needA
  if ~exist('r', 'var')
    d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]); % C-# displacements mat
    r = abs(d);                                    % dist matrix R^{MxN}
  end
  if self, r(diagind(r)) = 999; end % dummy nonzero diag values
end
if needA
  A = layerpot.fundsol(r, k);       % Phi
  Sker = A;                      % make kernel available w/out quadr w
else
  A = o.Sker; Sker = A;          % pass out Sker even if was passed in
end

if self % ........... source curve = target curve; can be singular kernel
  
  if isfield(o, 'self')
    A = o.self.S;                  % spit back the stored self-int mat
    
  elseif s.qtype=='p' & o.quad=='k'  % Kapur-Rokhlin (kills diagonal values)
    A(diagind(A)) = 0;
    [s w] = quadr.kapurtrap(N+1, o.ord);  % Zydrunas-supplied K-R weights
    w = 2*pi * w;                 % change interval from [0,1) to [0,2pi)
    A = circulant(w(1:end-1)) .* A .* repmat(sp.', [M 1]); % speed
    
  elseif s.qtype=='p' & o.quad=='m' % Martensen-Kussmaul (Kress MCM 1991)
    if isempty(s.kappa)
      error('cant do Martensen-Kussmaul quadr without s.kappa available!')
    end
    if k==0                         % Laplace SLP has ln sing
      S1 = -1/4/pi;                 % const M_1/2 of Kress w/out speed fac
      A = A - S1.*circulant(log(4*sin(pi*(0:N-1)/N).^2)); % A=D2=M_2/2 "
      A(diagind(A)) = -log(sp)/2/pi;        % diag vals propto curvature
    else
      if ~exist('r')
        d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]); % C-# displacements mat
        r = abs(d);                                    % dist matrix R^{MxN}
        r(diagind(r)) = 999; % dummy nonzero diag values
      end
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
    
  elseif s.qtype=='p' & o.quad=='a' % ---Alpert log-quadrature w/ rolling diag
    A = A .* repmat(s.w, [N 1]);  % use seg usual quadr weights away from diag
    A = quadr.alpertizeselfmatrix(A, k, s, @layerpot.Skernel, o); % note no 2pi factors
  
  else  % ------ self-interacts, but no special quadr, just use seg's
    % Use the crude approximation of kappa for diag, usual s.w off-diag...
    A = A .* repmat(s.w, [M 1]);  % use segment usual quadrature weights
    fprintf('warning: SLP crude self-quadr will be awful, low order\n')
    if k==0
      % ...
      A(diagind(A)) = 0;
    else
      eulergamma = -psi(1);  % put the below in for crude diag correction:
      A(diagind(A)) = s.w.*(1i/4-(eulergamma-1+log(k*s.w/4))/2/pi);
    end
  end
  
else % ............................ distant target curve, so smooth kernel
  
  A = A .* repmat(s.w, [M 1]);       % use segment quadrature weights
  
  % Now overwrite needed rows of A using very slow adaptive gauss quadrature...
  if isfield(o,'close') & k>0        % use adaptive quadr
    if ~isfield(o,'closeacc'), o.closeacc=1e-12; end
    rows = find(min(r,[],2)<o.close);   % which eval target pts need adaptive
    x = s.t;                            % source quadrature nodes in [0,1]
    w = zeros(size(x));
    for i=1:numel(x)       % Barycentric weights for Lagrange interp
      w(i) = 1./prod(x(find((1:numel(x))~=i))-x(i));
    end
    for i=rows'                      % must make a row vector for loop to work!
      if s.qtype=='g' | s.qtype=='c'
        for j=1:N
          xneqj = x(find((1:numel(x))~=j));  % col vec of nodes excluding j
          f = @(y) reshape(prod(repmat(xneqj, [1 numel(y)])-repmat(y(:).',[N-1 1]),1).' .* abs(s.Zp(y(:))) .* besselh(0,k*abs(s.Z(y(:))-t.x(i))), size(y));
          %f = @(y) y*fprintf('%g\n', y);  % debug tool to see what y requested
          %f([.2 .3])                      % debug
          A(i,j) = w(j) * (1i/4)*quadgk(f, 0, 1, 'RelTol', o.closeacc, 'AbsTol', o.closeacc, 'MaxIntervalCount', 1e3);
        end
      end
    end
  end
end
