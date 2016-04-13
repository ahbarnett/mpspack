function [t V] = gsvdtension(d, k, opts)
% GSVDTENSION  find various forms of tension for domain at given frequency
%
% t = gsvdtension(d, k, opts) returns tensions for domain d at wavenumber k.
%   t is a list of all the gen singular values, so min(t) is the tension.
%
% This is a development code for the ninc Neumann inclusion paper. It is much
%  cleaner and more accuracte than old hacks evp.tension, evp.tensionsq, and
%  finally gets the GSVD usage and eigenvectors correct.
%
% Options:
% opts.ten = 'v' : value part only of int Neu norm (t_simp), ie not true int nrm
%            'h' : Hassell suggestion use sqrt(G) for int nrm
%            'b' : Barnett realising can fill A,B first then regularize [A;B]
% opts.Fh = true, use Fh(1-h^2 \Delta_\pO) spectral weight op on u_n
% opts.solv: how linear algebra done: 'g' - GSVD, 'e' square then GEVP
%            (not implemented)
% opts.reg = 'n' : none
%            's' : SVD for new column space
%            't' : Timo's QR+SVD for new column space
% opts.eps = singular value relative cutoff for regularization
%
% [t V] = ... also returns N*N matrix of col sing vecs corresp to each singular
%   value.
%
% Notes: Makes obsolete @evp/tensionsq.m
% Allows various experimental versions of tension, for ninc B-Hassell-Tacy
% Contains helpers for numerical range of matrix.
% To do: 1) add in various Dirichlet tensions; 2) self-test, etc.
%
% See also: EVALBASES, EVP.INTNORMMATRIX

% Barnett 12/7/15 (my 43rd birthday)

if nargin<3, opts = []; end
if ~isfield(opts, 'ten'), opts.ten = 'v'; end
if ~isfield(opts, 'Fh'), opts.Fh = 0; end
if ~isfield(opts, 'solv'), opts.solv = 'g'; end
wantvecs = nargout>1;
ABon = 0;                                % by default, [A;B] not orthonormal.

E = k^2; d.k = k;
s = d.seg;   % the one segment
xn = real(conj(s.x) .* s.nx);        % x.n weight factor on bdry
w = s.w';    % bdry quadr weights as col vec

if opts.ten=='v'
  [A An] = d.evalbases(s);              % will get overwritten
  % regularize (transforms basis; need back-transform below)...
  [A V siginv N] = numrange([A;An],opts);  % note N (=rank) may be < Nf in d.bas
  [A An] = unstack(A,2);
  B = repmat(sqrt(xn .* w/2),[1 N]) .* A; % u-value part of Neu-case int norm
  if opts.Fh, An = evp.spectralfiltermatrix(s,k) * An; end   % apply bdry filter
  A = repmat(sqrt(w),[1 N]) .* An;        % Neu bdry norm
  
elseif opts.ten=='h'                      % Hassell idea to sqrt G
  [A Dx Dy] = d.evalbases(s);
  % regularize (transforms basis; need back-transform below)...
  [A V siginv N] = numrange([A;Dx;Dy],opts);  % N (=rank) may be < Nf in d.bas
  [A Dx Dy] = unstack(A,3);
  [G An] = evp.intnormmatrix(s, k, A, Dx, Dy);
  [W D] = eig(G);                         % get B as regularized sqrt of G:
  ev = diag(D); ii = (ev>=opts.eps*max(ev)); B = diag(sqrt(ev(ii)))*W(:,ii)';
  if opts.Fh, An = evp.spectralfiltermatrix(s,k) * An; end % apply bdry filter?
  A = repmat(sqrt(w),[1 N]) .* An;        % Neu bdry norm

elseif opts.ten=='b'                      % Barnett doing reg of [A;B] *after*
  [A Dx Dy] = d.evalbases(s);             %   .. building B=sqrt(G); careful
  [G An] = evp.intnormmatrix(s, k, A, Dx, Dy);
  [M N] = size(A);
  if opts.Fh, An = evp.spectralfiltermatrix(s,k) * An; end  % apply bdry filter?
  A = repmat(sqrt(w),[1 N]) .* An;        % Neu bdry norm "sqrt of F" matrix
  [W D] = eig(G);                         % get B as regularized sqrt of G
  ev = diag(D); ii = (ev>=opts.Geps*max(ev)); B = diag(sqrt(ev(ii)))*W(:,ii)';
  % regularize (transforms basis; need back-transform below)...
  [Q V siginv N] = numrange([A;B],opts);  % N (=rank) may be < Nf in d.bas
  A = Q(1:M,:); B = Q(M+1:end,:);         % note different heights!
  ABon = 1;                               % [A;B] = [Q_A;Q_B] are o.n.
end

if wantvecs                    % do GSVD
  [~,~,X,C,S] = gsvd(A,B);     % note X nonsingular, but not in general orthog
  t = sqrt(diag(C'*C)./diag(S'*S));
  if ~ABon                     % (if [A;B] o.n., X=X' from CS anyway, skip)
    X = inv(X');               % since Matlab's GSVD differs from Golub-vL!
    % (see gsvd help; annoying that have to do extra inversion step!)
  end
  V = V * (diag(siginv)*X);    % back xform gen sing vecs (since prefer 2-step)
else
  t = gsvd(A,B);               % faster if no vecs
end
%%%%%%%%%


function [Ar V siginv r] = numrange(A,opts)
% NUMRANGE - numerical range (column space) of matrix & its back-transformation
%
% [Ar V siginv r] = numrange(A,opts)
%
% Returns Ar : onb for numerical column space of A
%         V, siginv : parts of the back-transformation: if beta is vector that
%                  Ar would act on, then alpha = V * (diag(siginv)*beta) is the
%                  corresponding vector that A acts on. siginv is col vec.
%                  Note believe better passing 2 than combined matrix V*siginv.
%         r : rank (new # cols)
% opts.reg = 'n' : none
%            's' : SVD for new column space
%            't' : Timo's QR+SVD for new column space
% opts.eps = singular value relative cutoff for regularization
%
% Barnett 12/7/15
if ~isfield(opts, 'reg'), opts.reg = 't'; end
if ~isfield(opts, 'eps'), opts.eps = 1e-14; end
[M N] = size(A);
if strcmp(opts.reg,'n')
  Ar = A;
  V = eye(N); siginv = ones(N,1);
  r = N;
else                            % regularize (reg='s' or 't')
  if strcmp(opts.reg,'t'), [Q,A] = qr(A,0); end   % Timo's QR for speed (A<-R)
  [U S V] = svd(A, 'econ'); s = diag(S); clear S
  ii = (s >= opts.eps*max(s));  % boolean of cols to keep
  if strcmp(opts.reg,'t'), Ar = Q*U(:,ii); else, Ar = U(:,ii); end  % Timo
  V = V(:,ii); siginv = 1./s(ii);  % how to back transform vecs (finally ok!)
  r = sum(ii);                  % num rank
  %fprintf('  num rank = %d\n',r)
end

function [A B C] = unstack(Z,n)    % utility
if n>3, error('n at most 3'); end
[M N] = size(Z); m = round(M/n);
if m~=M/n, error('n must divide height of Z!'); end
A = Z(1:m,:);
if n>1, B = Z(m+1:2*m,:); end
if n>2, C = Z(2*m+1:3*m,:); end
