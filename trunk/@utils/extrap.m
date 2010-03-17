function [f0 err N] = extrap(f, hmax, opts)
% EXTRAP - Richardson extrapolate to evaluate f(0) via f(h), h>0, h dyadic mesh
%
% f0 = EXTRAP(f, hmax) where f is function handle, evaluates f(0) using f(h)
%   for h = hmax, hmax/2, hmax/4, ...  If too many function evals are required
%   a warning is given, suggesting hmax is descreased.
%
% f0 = EXTRAP(f, hmax, opts) allows options to be passed in:
%   opts.N - fixed number of evals (overrides accuracy considerations)
%   opts.reltol = desired relative accuracy (default 1e-14)
%   opts.abstol = desired absolute accuracy (default 1e-14)
%   opts.h0 = scale factor for rescaling polynomial basis (default 1)
%   opts.debug = true: print output of N, function values, current best guess
%  Note: the algorithm is satisfied when abstol *or* reltol is achieved.
%
% [f0 err] = ... also returns estimate of absolute error (difference of last
%   two cases as N is increased by 1)
%
% [f0 err N] = ... also returns the N (# function evals) used
%
% f may also return a row vector, in which case f0 and err are similarly-shaped.
% In this case, the worst-case of the elements of f is used for convergence.
%
% Notes: It seems to be good to choose hmax < (distance to singularity)/10.
% It's weird that polynomial basis does worse when scaled by hmax.
%
% Also see: UTILS/TESTEXTRAP for its test routine
%
% (C) Alex Barnett 3/15/10

warning('off', 'MATLAB:nearlySingularMatrix');
if nargin<3, opts = []; end
if ~isfield(opts,'reltol'), opts.reltol = 1e-14; end     % default rel acc
if ~isfield(opts,'abstol'), opts.abstol = 1e-14; end     % default abs acc
if ~isfield(opts,'h0'), opts.h0 = 1; end, h0=opts.h0;    % default h-scale
deb = isfield(opts, 'debug') && opts.debug;              % default is false

if isfield(opts,'N')               % fixed # evals
  N = opts.N;
  h = hmax*2.^-(0:N-1)';        % extrap pts mesh (col vec)
  fh = [];
  for j=1:N, fh = [fh; f(h(j))]; end        % data (col vec or matrix)
  V = zeros(N,N);
  for j=1:N, V(:,j) = (h/h0).^(j-1); end % col of V is poly basis eval on mesh
  c = V\fh;
  f0 = c(1,:);
  if deb
      fprintf('N=%d: h=%.3g\n', N, h(end));
      fprintf('\tf(h) = '), fprintf('%.16g ', fh(end,:));
      fprintf('\n\tf0 = '), fprintf('%.16g ', c(1,:)); fprintf('\n');
  end
  if nargout>1 % estimate error using difference from with one less mesh pt...
    h = h(1:end-1); fh = fh(1:end-1,:); V = V(1:end-1,1:end-1);
    c = V\fh; err = abs(f0-c(1,:));
  end

else                   % desired rel err; steadily increase N until good enough
  Nmax = 14; f0 = nan; h = []; fh = [];   % Max N here; build up array as go
  for N=1:Nmax
    h = [h; hmax*2^-(N-1)]; fh = [fh; f(h(end))]; M = size(fh, 2);
    V = zeros(N,N);    % Vandermonde type matrix; note poly basis not scaled!
    for j=1:N, V(:,j) = (h(1:N)/h0).^(j-1); end % col of V = poly eval on mesh
    c = V\fh;
    if deb
      fprintf('N=%d: h=%.3g\n', N, h(end));
      fprintf('\tf(h) = '), fprintf('%.16g ', fh(end,:));
      fprintf('\n\tf0 = '), fprintf('%.16g ', c(1,:)); fprintf('\n');
    end
    err = abs(c(1,:)-f0);
    if max(err./abs(c(1,:)))<=opts.reltol || max(err)<=opts.abstol, ...
          f0 = c(1,:); break; end
    f0 = c(1,:);
  end
  if N==Nmax,warning 'extrap hit max # func evals (20); try reducing hmax';end
end
warning('on', 'MATLAB:nearlySingularMatrix');  % or return to global state?
