function [A Dker_noang cosker] = D(k, s, t, o)
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
%                'a' (Alpert)
%                periodic quadrature rules, used only if s.qtype is 'p';
%                any other does low-order non-periodic quad using segment's own.
%    opts.ord = 2,6,10,etc controls order of Kapur-Rokhlin or Alpert rule.
%    opts.Dker_noang = quad-unweighted kernel matrix of fund-sols derivs
%                      without the cosphi factors (prevents recomputation
%                      of H1's).
%    opts.displ = matrix of source-target displacements, prevents recalculation.
%    opts.rdist = matrix of source-target distances, prevents recalculation.
%    opts.derivSLP = true, then uses target-normals instead of
%                      source-normals, which computes n-deriv of SLP instead.
%    opts.close = distance below which adaptive quadrature is used to evaluate
%                distant targets (slow). If not present, never adaptive.
%    opts.closeacc = set relative tolerance for close eval (default 1e-12)
%
%  [D Dker_noang cosker] = D(...) also returns quad-unweighted kernel values
%    matrix Dker_noang, and matrix of cos angle factors (cosphi or costh)
%
% Issues: need to remove the diag=999 hack (affects slow matlab eval only)!

% Copyright (C) 2008 - 2012, Alex Barnett and Timo Betcke

if isempty(k) | isnan(k), error('DLP: k must be a number'); end
if nargin<4, o = []; end
if nargin<3, t = []; end
if ~isfield(o, 'quad'), o.quad='m'; end; % default self periodic quadr
if ~isfield(o, 'ord'), o.ord=10; end;    % default quadr order
symmflagval = -999;   % all diag vals of this signifies symmetric - a hack
self = isempty(t);               % self-interact: potentially sing kernel
N = numel(s.x);                  % # src pts
if self, t = s; end              % use source as target pts (handle copy)
M = numel(t.x);                  % # target pts
sp = s.speed/2/pi; % note: 2pi converts to speed wrt s in [0,2pi] @ src pt
needA = ~isfield(o, 'Dker_noang'); % true if must compute H1 matrix
if isfield(o, 'displ'), d = o.displ; else
  d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]); end % C-# displacements mat
if isfield(o, 'rdist'), r = o.rdist; else
  r = abs(d); end                                    % dist matrix R^{MxN}
if self, r(diagind(r)) = symmflagval; end % flag diag vals as self-int
dSLP=0; if isfield(o, 'derivSLP') & o.derivSLP, dSLP=1; end  % flag
if dSLP                           % dPhi(x,y)/dnx, ie D^T adjoint op
  nx = repmat(-t.nx, [1 N]);      % identical cols given by targ normals
else                              % dPhi(x,y)/dny, ie D op
  nx = repmat(s.nx.', [M 1]);     % identical rows given by src normals
end
cosker = real(conj(nx).*d) ./ r;  % dot prod <normal, displacement>
if needA
  [A Dker_noang] = layerpot.fundsol_deriv(r, cosker, k); % n-deriv of Phi
else
  [A Dker_noang] = layerpot.fundsol_deriv(r, cosker, k, o.Dker_noang);
end                               % A is now the kernel value matrix

if self % ........... source curve = target curve; can be singular kernel

  if s.qtype=='p' & o.quad=='k' % Kapur-Rokhlin (kills diagonal values)
    A(diagind(A)) = 0;
    [s w] = quadr.kapurtrap(N+1, o.ord);  % Zydrunas-supplied K-R weights
    w = 2*pi * w;                 % change interval from [0,1) to [0,2pi)
    A = circulant(w(1:end-1)) .* A .* repmat(sp.', [M 1]); % speed
    
  elseif s.qtype=='p' & o.quad=='m' % Martensen-Kussmaul (Kress MCM 1991)
    if isempty(s.kappa)
      error('cant do Martensen-Kussmaul quadr without s.kappa available!')
    end
    if k==0                   % for k=0 DLP analytic on analytic curve
      A(diagind(A)) = -s.kappa/(4*pi);        % diag vals propto curvature
      A = (2*pi/N) * A .* repmat(sp.', [M 1]);   % speed quad weights
    else                      % for k>0 has x^2 ln x sing, Kress handles
                              % (note, diag A is garbage right now)
      D1 = triu(besselj(1,k*triu(r,1)),1); % use symmetry (arg=0 is fast)
      D1 = -(k/4/pi)*cosker.*(D1.'+D1);  % L_1/2 of Kress w/out speed fac
      A = A - D1.*circulant(log(4*sin(pi*(0:N-1)/N).^2)); % A=D2=L_2/2 "
      A(diagind(A)) = -s.kappa/(4*pi);   % L_2(t,t)/2, same as for k=0
      % speed factors: diag matrix mult from right...
      A = (circulant(quadr.kress_Rjn(N/2)).*D1 + (2*pi/N)*A) .* ...
          repmat(sp.', [M 1]);
    end
    
  elseif s.qtype=='p' & o.quad=='a' % ---Alpert log-quadrature w/ rolling diag
    A = A .* repmat(s.w, [N 1]);  % use seg usual quadr weights away from diag
    kerfun = @layerpot.Dkernel; if dSLP, kerfun = @layerpot.DTkernel; end
    A = quadr.alpertizeselfmatrix(A, k, s, kerfun, o);  % note no 2pi factors
  
  else       % ------ self-interacts, but no special quadr, just use seg's
    % Use the crude approximation of kappa for diag, usual s.w off-diag...
    if isempty(s.kappa)
      fprintf('warning: DLP crude self-quadr has no s.kappa so will be awful!\n')
      A(diagind(A)) = 0;
    else
      if k~=0, fprintf('warning: DLP crude self-quadr, using s.kappa for diag, only good for k=0\n'); end
      A(diagind(A)) = -s.kappa/(4*pi);  % diag vals propto curvature
      % (correct for k=0 Laplace only)
    end
    A = A .* repmat(s.w, [M 1]);  % use segment usual quadrature weights
  end
  
else % ............................ distant target curve, so smooth kernel
  
  A = A .* repmat(s.w, [M 1]);       % use segment quadrature weights

  % Now overwrite needed rows of A using very slow adaptive gauss quadrature...
  if 1 & isfield(o,'close') & k>0        % use adaptive quadr
    if ~isfield(o,'closeacc'), o.closeacc=1e-12; end
    rows = find(min(r,[],2)<o.close);   % which eval target pts need adaptive
   %numel(rows)
    x = s.t;                            % source quadrature nodes in [0,1]
    w = zeros(size(x));
    for i=1:numel(x)       % Barycentric weights for Lagrange interp
      w(i) = 1./prod(x(find((1:numel(x))~=i))-x(i));
    end
    for i=rows'                      % must make a row vector for loop to work!
      if s.qtype=='g' | s.qtype=='c'
        for j=1:N                    % loop over Lagrange basis funcs for src
          xneqj = x(find((1:numel(x))~=j));  % col vec of nodes excluding j
          if dSLP
            f = @(y) Lagrange_SLP_deriv(y, xneqj, k, s, t.x(i), t.nx(i));
          else
            f = @(y) Lagrange_DLP(y, xneqj, k, s, t.x(i));
          end
%          fprintf('basis func %d:\n', j)
          A(i,j) = w(j) * (1i*k/4)*quadgk(f, 0, 1, 'RelTol', o.closeacc, 'AbsTol', o.closeacc, 'MaxIntervalCount', 1e3);
        end
      end
    end
  end
end

function f = Lagrange_SLP_deriv(t, xneqj, k, s, x, nx) % ---------------------
% evaluate target-deriv of j-th Lagrange basis density SLP
% at any list of t (0<=t<=1 along source segment).
% x, nx = single target location and normal deriv.
%
% Used by slow adaptive close eval
% Lagrange basis (excluding its t-indep denominators)...
f = prod(repmat(xneqj, [1 numel(t)])-repmat(t(:).',[numel(xneqj) 1]), 1).';
d = s.Z(t(:))-x; r = abs(d);
% now mult by SLP deriv func...
f = f .* abs(s.Zp(t(:))) .* besselh(1,k*r);
f = f .* real(conj(nx) * d) ./ r;            % ? use target -nx (cf start above)
f = reshape(f, size(t));


function f = Lagrange_DLP(t, xneqj, k, s, x) % -------------------------------
% evaluate value of j-th Lagrange basis density DLP
% at any list of t (0<=t<=1 along source segment).
% x = single target location. Used by slow adaptive close eval.
%
% Lagrange basis (excluding its t-indep denominators)...
f = prod(repmat(xneqj, [1 numel(t)])-repmat(t(:).',[numel(xneqj) 1]), 1).';
d = s.Z(t(:))-x; r = abs(d);
% now mult by SLP deriv func...
f = f .* abs(s.Zp(t(:))) .* besselh(1,k*r);
f = f .* real(conj(-s.Zn(t(:))) .* d) ./ r;  % opposite sign for src normals
f = reshape(f, size(t));
%fprintf('\tquadgk using %d pts\n', numel(t))
