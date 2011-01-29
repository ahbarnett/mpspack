function A = alpertizeselfmatrix(A, k, s, kerfun, o) 
% A = alpertizeselfmatrix(A, k, s, kerfun, opts) returns Alpert quadrature rule
%  matrix corrected along a band around the diagonal, given the matrix A of
%  source-to-target kernel values (including speed function and uniform
%  quadrature weights 2pi/N). k is omega the wavenumber, kerfun is a
%  kernel function of the form val = kerfun(k, x, nx, y, ny)
%  returning k(s,t) without the speed |z'(t)| factor.
%  s is the segment whose self-interaction is needed.
%  opts.ord controls Alpert quadrature order; other options ignored.
%  Now vectorized for speed! (Most time is spent on kernel evaluation, good)
%
% Issues: Matlab makes local copy of whole of A, which is memory-movement bad!
%
% Also see: LAYERPOT.S, LAYERPOT.D, LAYERPOT.T, TEST/TESTLPQUAD

% 1/29/11. Copyright (C) 2011 Alex Barnett
M = size(A,1); N = size(A,2);
if M~=N, warning('Alpert quadr will fail unless matrix square!'); end
[tex,wex,nskip] = quadr.QuadLogExtraPtNodes(o.ord); % get log nodes, wei
tex = [-tex(end:-1:1); tex]; wex = [wex(end:-1:1); wex]; % make symmetric
wex = wex/N;      % scale for curve regular node parameter spacing
Na = numel(tex);   % # Alpert pts (both sides)
if N<=2*nskip, warning('not enough regular nodes for Alpert order!'); end
ninterp = o.ord + 3;   % # local regular Lagrange interpolation pts (O'Neil)

for i=1:N     % loop over rows (i=targets, j=source density dofs)
  A(i, mod(i+[-nskip+1:nskip-1]-1, N)+1) = 0; % kill skipped j's near diag
end

% These arrays are N-by-Na...
x = repmat(s.x, [1 Na]); nx = repmat(s.nx, [1 Na]); % target locations & normals
t = repmat(s.t, [1 Na]) + repmat(tex.'/N, [N 1]); % block src curve parameters
% Get all curve info from segment at once = fast! :
y = s.Z(t); ny = s.Zn(t); sp = abs(s.Zp(t));  % speed funcs at src pts
K = sp .* kerfun(k, x, nx, y, ny); % get all kernel evals at once = fast
% Note chould change the above for quasi-periodic open segments...

for l=1:numel(tex)  % loop over Alpert's nonregular quadr pts
  dt = tex(l)/N;    % parameter offset for current quadr pt
  joffs = floor(tex(l)-ninterp/2+1) + (0:ninterp-1); % nearest interp pts
  % interpolation pts and L = row vec of Lagrange basis funcs eval at dt:
  x = joffs/N; w = utils.baryweights(x); L = utils.baryprojs(x, w, dt);
  ii = kron(ones(1,ninterp), 1:N);  % i-indices of A to write over
  jj = mod(ii + kron(joffs, ones(1,N))-1, N) + 1;  % j-indices ditto
  iii = sub2ind(size(A), ii,jj);      % 1d array indices ditto
  A(iii) = A(iii) + wex(l) * kron(ones(1,ninterp), K(:,l).') .* kron(L, ones(1,N));
end
% Note the affected elements spill out to dist nskip-1+iord/2 from diag!
% (even though the diag band that was set to zero is only to nskip dist)
