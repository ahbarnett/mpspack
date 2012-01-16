function [A Sker Dker_noang] = T(k, s, t, o);
% T - deriv of double layer potential discr matrix for density on a segment
%
%  T = T(k, s) where s is a segment returns the matrix discretization of
%   the deriv-of-DLP integral operator with density on the segment,
%      u(x) = int_s partial_n(x) partial_n(y) Phi(x,y) tau(y) ds(y)
%   where Phi is the fundamental solution
%      Phi(x,y) = (i/4) H_0^{(1)}(k|x-y|)    for k>0
%                 (1/2pi) log(1/|x-y|)       for k=0
%   for x each of the points on the self-same segment s. No jump relations
%   are included. For self-interaction, no special quadrature is yet implemented
%
%  T = T(k, s, t) where t is a pointset (any object with fields t.x)
%   uses s as the source segment as above, but target points in t.
%   It is assumed t is not s.
%
%  T = T(k, s, [], opts) or T = T(k, s, t, opts) does the above two choices
%   but with options struct opts including the following:
%    opts.quad = 'k' (Kapur-Rokhlin), 'm' (Maue-Kress spectral* for k>0)
%                'a' (Alpert)
%                periodic quadrature rules, used only if s.qtype is 'p';
%                any other does low-order non-periodic quad using segment's own.
%      * Note: the k-independent cot (Cauchy pole) is dropped since cancels
%        (T-T_0) in Rokhlin scheme [Kress 1991 MCM gives needed Wittich quadr].
%    opts.ord = 2,6,10,etc controls order of Kapur-Rokhlin or Alpert rule.
%    opts.Sker = quad-unweighted kernel matrix of fund-sols (i/4)H_0
%    opts.Dker_noang = quad-unweighted kernel matrix of fund-sols derivs
%                      without the cosphi factors (prevents recomputation
%                      of H1's), ie (ik/4)H_1
%    opts.displ = matrix of source-target displacements, prevents recalculation.
%    opts.rdist = matrix of source-target distances, prevents recalculation.
%   Passing any of the above options stops it from being calculated internally.
%   When k=0, the options Sker and Dker_noang have no effect.
%    opts.close = distance below which adaptive quadrature is used to evaluate
%                distant targets (slow). If not present, never adaptive.
%    opts.closeacc = set relative tolerance for close eval (default 1e-12)
%
%  [T Sker Dker_noang] = T(...) also returns quad-unweighted
%   kernel values matrices Sker, Dker_noang, when k>0 (empty for Laplace)
%

% Copyright (C) 2008 - 2012, Alex Barnett and Timo Betcke


if isempty(k) | isnan(k), error('T: k must be a number'); end
if nargin<4, o = []; end
if nargin<3, t = []; end
if ~isfield(o, 'quad'), o.quad='m'; end; % default self periodic quadr
if ~isfield(o, 'ord'), o.ord=10; end;    % default quadr order (Kapur-Rokh)
symmflagval = -999;   % all diag vals of this signifies symmetric - a hack
self = isempty(t);               % self-interact: potentially sing kernel
N = numel(s.x);                  % # src pts
if self, t = s; end              % use source as target pts (handle copy)
M = numel(t.x);                  % # target pts
sp = s.speed/2/pi; % note: 2pi converts to speed wrt s in [0,2pi] @ src pt
% compute all geometric parts of matrix...
if isfield(o, 'displ'), d = o.displ; else
  d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]); end % C-# displacements mat
if isfield(o, 'rdist'), r = o.rdist; else
  r = abs(d); end                                    % dist matrix R^{MxN}
if self, r(diagind(r)) = symmflagval; end % dummy diag values flag self
ny = repmat(s.nx.', [M 1]);      % identical rows given by src normals
csry = conj(ny).*d;              % (cos phi + i sin phi).r
nx = repmat(t.nx, [1 N]);        % identical cols given by target normals
csrx = conj(nx).*d;              % (cos th + i sin th).r
clear d
if ~isfield(o, 'Sker')
  o.Sker = layerpot.fundsol(r, k);           % Phi
end
if ~self | s.qtype~='p' | o.quad~='m' % all but hypersing spectral quadr
  clear nx ny
  if k==0                          % Laplace; don't reuse since v. fast anyway
    A = -real(csry.*csrx)./(r.^2.^2)/(2*pi); % (-1/2pi)cos(phi-th)/r^2
    Sker = []; Dker_noang = [];
  else
    cc = real(csry).*real(csrx) ./ (r.*r);      % cos phi cos th
    cdor = real(csry.*csrx) ./ (r.*r.*r);   % cos(phi-th) / r
    clear csrx csry
    if ~isfield(o, 'Dker_noang')
      [A o.Dker_noang] = layerpot.fundsol_deriv(r, -cdor, k); % get n-deriv Phi
    else
      [A o.Dker_noang] = layerpot.fundsol_deriv(r, -cdor, k, o.Dker_noang);
    end
    % put all this together to compute the kernel value matrix...
    A = A + k^2 * cc .* o.Sker;
    Sker = o.Sker; Dker_noang = o.Dker_noang; % pass out some calcs for reuse
  end
end

if self % ........... source curve = target curve; can be singular kernel
  
  if isfield(o, 'self')
    A = o.self.T;                  % spit back the stored self-int mat
    
  elseif s.qtype=='p' & o.quad=='k' % Kapur-Rokhlin (ignores diagonal value)
    A(diagind(A)) = 0;
    [s w] = quadr.kapurtrap(N+1, o.ord);  % Zydrunas-supplied K-R weights
    w = 2*pi * w;                 % change interval from [0,1) to [0,2pi)
    A = circulant(w(1:end-1)) .* A .* repmat(sp.', [M 1]); % speed
  
  elseif s.qtype=='p' & o.quad=='m' % hypersing spectral, Maue '49, Kress '91
    if isempty(s.kappa) | isempty(s.relaccel)
      error('cant do Maue-Kress hypersing quadr without s.kappa or s.relaccel!')
    end
    if k==0
      error 'Maue-Kress hypersing quadr not implemented for k=0!'
    else
      if ~exist('r', 'var')  % compute dist matrix if needed
        d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]); % C-# displacements mat
        r = abs(d);                                    % dist matrix R^{MxN}
        r(diagind(r)) = 999; % dummy nonzero diag values
      end
      AS = o.Sker;                       % single-layer kernel value matrix
      S1 = triu(besselj(0,k*triu(r,1)),1);  % use symmetry (arg=0 is fast)
      S1 = -(1/4/pi)*(S1.'+S1);     % next fix it as if diag(r) were 0
      S1(diagind(S1)) = -(1/4/pi);  % S1=M_1/2 of Kress w/out speed fac
      logmat = circulant(log(4*sin(pi*(0:N-1)/N).^2));
      AS = AS - S1.*logmat;                   % AS=D2=M_2/2 "
      eulergamma = -psi(1);         % now set diag vals Kress M_2(t,t)/2
      AS(diagind(AS)) = 1i/4 - eulergamma/2/pi - log((k*sp).^2/4)/4/pi;
      AS = (circulant(quadr.kress_Rjn(N/2)).*S1 + AS*(2*pi/N)) .* ...
          repmat(sp.', [M 1]);
      % ...AS is now single-layer Kress-quadr matrix w/ speed factor
      A = k^2 * real(conj(nx).*ny) .* AS;   % first Maue term
      % Do tangential deriv term (2nd term in Maue splitting eqn):
      spsinphi = imag(csrx)./r .* repmat(sp, [1 M]);  % -speed times sin phi 
      if ~isfield(o, 'Dker_noang')
        [B o.Dker_noang] = layerpot.fundsol_deriv(r, -spsinphi, k);
      else
        [B o.Dker_noang] = layerpot.fundsol_deriv(r, -spsinphi, k,o.Dker_noang);
      end
      S1 = triu(besselj(1,k*triu(r,1)),1);  % use symmetry (arg=0 is fast)
      S1 = (k/4/pi)*(S1.'+S1) .* spsinphi;  % S1=N_1/2 of Kress, diag is 0
      B = B - circulant(cot(pi*(0:N-1)/N)/4/pi) - S1.*logmat; % B=N_2/2 Kress
      B(diagind(B)) = -s.relaccel/4/pi/2/pi; % diag val Kress, sign err? (/2pi) 
      %figure; imagesc(real(B)); colorbar; figure; plot(real(B(70,:)), '+-');
      D = circulant(quadr.perispecdiffrow(N)); % spectral differentiation
      % add 2nd term in Maue...
      A = A + (repmat(1./sp, [1 M]).*(circulant(quadr.kress_Rjn(N/2)).*S1 + B*(2*pi/N))) * D;
    end
    
  elseif s.qtype=='p' & o.quad=='a' % ---Alpert log-quadrature w/ rolling diag
    A = A .* repmat(s.w, [N 1]);  % use seg usual quadr weights away from diag
    A = quadr.alpertizeselfmatrix(A, k, s, @layerpot.Tkernel, o); % note no 2pi factors
    % NOTE: Alpert only applies log-sing correctly, not the hypersingular T,
    % but if two T's are subtracted, it will work fine for the difference.
  
  else       % ------ self-interacts, but no special quadr, just use seg's
    % Use the crude approximation of kappa for diag, usual s.w off-diag...
    A = A .* repmat(s.w, [M 1]);  % use segment usual quadrature weights
    warning 'using T crude self-quadr, will be awful (no diag)!'
    A(diagind(A)) = 0;
  end
  
else % ............................ distant target curve, so smooth kernel
  
  A = A .* repmat(s.w, [M 1]);       % use segment quadrature weights
  
  % Now overwrite `close' rows of A, if any, using very slow adaptive gauss quadr...
  if isfield(o,'close') & k>0        % use adaptive quadr
    if ~isfield(o,'closeacc'), o.closeacc=1e-12; end
 %r   
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
          f = @(y) Lagrange_DLP_deriv(y, xneqj, k, s, t.x(i), t.nx(i));
          A(i,j) = w(j) * quadgk(f, 0, 1, 'RelTol', o.closeacc, 'AbsTol', o.closeacc, 'MaxIntervalCount', 1e4);
        end
      end
    end
  end
end

function f = Lagrange_DLP_deriv(t, xneqj, k, s, x, nx) % -----------------------
% evaluate target-normal deriv of j-th Lagrange basis density DLP
% at any list of t (0<=t<=1 along source segment).
% x, nx = single target location, normal deriv.
%
% Lagrange basis (excluding its t-indep denominators)...
f = prod(repmat(xneqj, [1 numel(t)])-repmat(t(:).',[numel(xneqj) 1]), 1).';
d = s.Z(t(:))-x; r = abs(d);
csrx = conj(nx)*d;                     % (code taken from above)
csry = conj(s.Zn(t(:))).*d;             % cos src normals
cc = real(csry).*real(csrx) ./ (r.*r);      % cos phi cos th
cdor = real(csry.*csrx) ./ (r.*r.*r);   % cos(phi-th) / r
val = (1i*k/4)*besselh(1,k*r) .* (-cdor) + (1i*k*k/4)*cc.*besselh(0,k*r);
f = reshape(f .* val .* abs(s.Zp(t(:))), size(t));
%fprintf('%.16g\n', t);                 % for diagnostics of quadgk