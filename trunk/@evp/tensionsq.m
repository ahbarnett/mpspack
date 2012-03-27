function [t V F G] = tensionsq(d, E, opts)
% TENSIONSQ - return tension eigen/singular values for domain at given frequency
%
% t = tensionsq(d, E, opts) finds tension squared for domain d,
%   at Helmholtz parameter (frequency squared) E=omega^2.
%   min(t) is tension squared; t are generalized eigenvalues of (F,G).
%
% [t V] = ... returns eigenvectors too as in Matlab's eig.
%
% opts struct may contain:
% meth = 'g' ('gg', 'gt', 'gs') : Timo's GSVD w/ Dirichlet approx interior norm.
% meth = 'i' uses interior points, needs:
%          opts.i = interior pointset, opts.a = area assoc to each int pt
% meth = 'a' ('at' etc) : correct interior norm via Alex formula, GEVP.
%
% gep : if GEVP used, controls way eig(...) is done.
%
% neu : if true, switch from Dirichlet to Neumann BCs (neu=1 basic, =2 F_h vers)
%       needs opts.arcl, row vec approx to cumulative arclength samples
%
% Notes: 1) This is a helper routine for MPS eigenvalue papers, particularly
% the PSPM conference paper with Andrew Hassell.
% 2) d must contain one segment, one basis, currently.
% 3) Very alpha code.
%
% A. Barnett, April 2010.

if nargin<3, opts = []; end
if ~isfield(opts, 'meth'), opts.meth = 'a'; end
if ~isfield(opts, 'gep'), opts.gep = 'v'; end
if ~isfield(opts, 'tensstar'), opts.tensstar = 0; end
if ~isfield(opts, 'neu'), neu=0; else neu=opts.neu; end % default neu=0
wantvecs = nargout>1;

s = d.seg;                           % get the one segment
N = d.bas{1}.Nf;
xn = real(conj(s.x) .* s.nx);        % x.n weight factor on bdry
d.k = sqrt(E);
wei = 1+0*xn; if opts.tensstar, wei = 1./xn; end % weighting function for A, F

if opts.meth(1)=='g'| opts.meth(1)=='i'   % GSVD: Timo's. overrides gep method
  if opts.meth(1)=='i'          % interior pts
    A = repmat(sqrt(s.w.'),[1 N]) .* d.evalbases(s);
    B = sqrt(opts.a) * d.evalbases(opts.i);
  else                      % Dirichlet-interior norm approx
    [A An] = d.evalbases(s);
    A = repmat(sqrt(s.w.'.*wei),[1 N]) .* A;           % bdry norm
    B = repmat(sqrt(xn .* s.w.' / (2*E)),[1 N]) .* An; % n-deriv approx int norm
  end
  if opts.meth(2)=='g'                 % don't regularize
    if wantvecs
      [UU VV V C S] = gsvd(A,B);  % V gives eigvecs
      t = sqrt(diag(C'*C)./diag(S'*S));
    else, t = gsvd(A,B); end
    
  elseif opts.meth(2)=='t'    % Timo's code; uses QR since faster than tall SVD
    [Q,R]=qr([A;B],0); [U,S,V]=svd(R); S=diag(S);
    ii = abs(S)>opts.eps*max(abs(S));   % note I scaled it to max(sig val)
    Q=Q*U(:,ii);  % cols of Q now onb for Col[A;B] at numerical rank.
    if wantvecs
      [UU VV X C S] = gsvd(Q(1:size(A,1),:),Q(size(A,1)+1:end,:));
      V = V(:,ii) * X;          % rotate eigvecs back (no R needed! weird)
      t = sqrt(diag(C'*C)./diag(S'*S));
    else, t = gsvd(Q(1:size(A,1),:),Q(size(A,1)+1:end,:)); end
    
  elseif opts.meth(2)=='s'      % use SVD to get numerical col space, slower
    [Q S V] = svd([A;B], 'econ'); S=diag(S);
    ii = abs(S)>opts.eps*max(abs(S)); Q=Q(:,ii);
    if wantvecs
      [UU VV X C S] = gsvd(Q(1:size(A,1),:),Q(size(A,1)+1:end,:));
      V = V(:,ii) * X;             % rotate eigvecs back
      t = sqrt(diag(C'*C)./diag(S'*S));
    else, t = gsvd(Q(1:size(A,1),:),Q(size(A,1)+1:end,:)); end
  end
  t = t.^2;  % square the gsingvals to match the GEP version
  return
  
elseif opts.meth(1)=='r'              % Dirichlet-approx to int norm (rellich)
  [A An] = d.evalbases(s);
  A = repmat(sqrt(s.w.'),[1 N]) .* A;           % bdry norm
  B = repmat(sqrt(xn .* s.w.' / (2*E)),[1 N]) .* An; % n-deriv approx int norm
  if opts.meth(2)=='t'
    [Q,R]=qr([A;B],0); [U,S,V]=svd(R); S=diag(S); % 
    ii = abs(S)>opts.eps*max(abs(S));   % note I scaled it to max(sig val)
    Q=Q*U(:,ii);  % cols of Q now onb for Col[A;B] at numerical rank.
    A = Q(1:size(A,1),:); B = Q(size(A,1)+1:end,:);
  end
  F = A'*A; G = B'*B;
  %F = (F+F')/2; G = (G+G')/2;   % no difference
  
elseif opts.meth(1)=='a'     % full interior norm computed via Alex formula
  [A Dx Dy] = d.evalbases(s);
  if length(opts.meth)>1 & opts.meth(2)=='t'  % Timo-style regularize via orthog
    [Q,R]=qr([A;Dx;Dy],0); [U,S,V]=svd(R); S=diag(S); % onb from all matrices
    ii = abs(S)>opts.eps*max(abs(S));   % note I scaled it to max(sig val)
    Q=Q*U(:,ii); N=size(Q,2); % cols of Q onb for Col[A;Ax;Ay] at num rank.
    M=size(A,1); A = Q(1:M,:); Dx = Q(M+1:2*M,:); Dy = Q(2*M+1:end,:);
  end
  Dnor = repmat(real(s.nx), [1 N]).*Dx + repmat(imag(s.nx), [1 N]).*Dy;
  if neu                         % Neumann BC case, new
    M=numel(s.w); L = sum(s.w); % perim. outerprod for Fourier proj matrix...
    ns = -M/2+1:M/2; P = exp((-2i*pi/L)*ns'*opts.arcl); % ns=freqs
    Pt = P'; P = P.*repmat(s.w/L,[M 1]); % Pt takes F coeffs to vals on pO
    %sum(P,2)  % should be zero apart from mode m=0 row which is 1
    %figure; semilogy(abs(P*exp(0*1i*real(s.x)))); % test exp decay of F coeffs
    if neu==1, F = Dnor'*(repmat(s.w.',[1 N]).*Dnor); % naive neumann
    else, k=sqrt(E); xin = 2*pi/L/k*ns; % overall wavenumber, xi values
      Fk=max(real(sqrt(1-xin.^2)),k^(-1/3)); % figure; plot(xin, Fk, '+-');
      F = diag(1./Fk)*P*Dnor; F = F'*F; % l^2 inner prod in Fourier coeffs
    end
  else   % usual Dirichlet
    F = A'*(repmat(s.w.',[1 N]).*A); % build F, then G adding pieces together...
  end
  W = repmat(s.w.'.*xn, [1 N]);              % x.n as elementwise mult to right
  G = A' * (W .* A) * E; clear A
  G = G - Dx' * (W .* Dx) - Dy' * (W .* Dy);
  W = repmat(s.w.', [1 N]);                  % 1 as elementwise mult to right
  Ddil = repmat(real(s.x), [1 N]).*Dx + repmat(imag(s.x), [1 N]).*Dy;
  B = Ddil'*(W.*Dnor); clear Ddil Dnor W
  G = G + B + B';                         % conjugate saves another mat mult
  G = (G+G')/(4*E);                       % rescale and make symm
end
  
if wantvecs
  if opts.gep=='e'                          % std
    [X D] = eig(F,G); t = diag(D);
  elseif opts.gep=='c'                      % cholesky, ok for v small systems
    [X D] = eig(F,G, 'qz'); t = diag(D);
  elseif opts.gep=='v'                      % vergini regularization method
    [t X] = utils.regeig(F,G, opts);
  elseif opts.gep=='i'                      % reversed, reg method
    [t X] = utils.regeig(G,F, opts); t = 1./t;
  end
  if opts.meth(1)~='g' & length(opts.meth)>1 & (opts.meth(2)=='t' | opts.meth(2)=='s')
    V = V(:,ii) * diag(1./S(ii)) * X;        % un-orthogonalize the eigvecs
  else
    V = X;
  end

else
  if opts.gep=='e'                          % std
    t = eig(F,G);
  elseif opts.gep=='c'                      % cholesky, ok for v small systems
    t = eig(F,G, 'qz');
  elseif opts.gep=='v'                      % vergini regularization method
    t = utils.regeig(F,G, opts);
  elseif opts.gep=='i'                      % reversed, reg method
  t = 1./utils.regeig(G,F, opts);
  end
end
