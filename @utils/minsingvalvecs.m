function [u s v info] = minsingvalvecs(A, opts)
% MINSINGVALVECS - iterative estimate of minimum singular value of square matrix
%
% [u s v info] = MINSINGVALVECS(A) returns left- and right-singular vectors u,v
%   and corresponding singular value s which is the minimum one for A.
%   info.flag = 0 indicates success, 1 failure. info.its = # its used (large
%   indicates s_{N-1} and s_{N} are close to each other, and results inaccurate.
%   For well-separated singular values, <10 its is enough, <3 if v. separated).
%
% [u s v info] = MINSINGVALVECS(A, opts) allows control of method parameters:
%    opts.maxits : maximum number of iterations (default 100)
%    opts.tol : acceptable relative tolerance in singular value (default 1e-14)
%
% Uses simplification of: Loef algorithm as outlined in Calderon-Guizar et al,
% IEEE Power Engineering Review, 19 (9), 55-56 (1999). This involves LU-decomp
% of A which is O(N^3) but overall up to 30 times faster than dense complex svd!
%
% Notes/issues:
%  * When s is small (close to 1e-16) the phase angle of u relative to v
%    becomes inaccurate. This may be an intrinsic property of SVD?
%  * Explore rectangular A variant (do LU of A' also, doubling the cost?)
%  * Make a block version to handle subspace of several small sing vals?
%
% See also: SVD, TEST/TESTMINSINGVALVEC

% Copyright (C) 2010, Alex Barnett
  if nargin<2, opts = []; end
  if ~isfield(opts, 'maxits'), opts.maxits = 100; end
  if ~isfield(opts, 'tol'), opts.tol = 1e-14; end 
  N = size(A,1); if N~=size(A,2), error('A must be square!'); end
 
  [L U] = lu(A); LT = L'; UT = U';   % the slow part; replace with permuted?
  v = rand(N,1) - 0.5;               % starting vector (real if A is)
  if ~isreal(A), v = v + 1i*(rand(N,1) + 0.5); end  % complex starting choice
  s = nan;                           % dummy starting for min sing val
  for i=1:opts.maxits, os = s;       % store the old value of s
    z=UT\v; u=LT\z; u=u/norm(u); z=L\u; v=U\z; s=1/norm(v); v=s*v;
    if abs((s-os)/os)<opts.tol, break; end     % stop if small relative change
  end
  r = abs(v(1))/v(1); v=v*r; u=u*r;  % make first el real, as matlab svd
  info.flag = (i==opts.maxits);      % error handling
  info.its = i;
