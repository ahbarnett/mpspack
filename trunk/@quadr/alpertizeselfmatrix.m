function A = alpertizeselfmatrix(A, k, s, kerfun, o) 
% A = alpertizeselfmatrix(A, k, s, kerfun, opts) returns Alpert quadrature rule
%  matrix corrected along a band around the diagonal, given the matrix A of all
%  source-to-target kernel values (including speed function and uniform
%  quadrature weights 2pi/N), with diag(A)=0. k is omega wavenumber, kerfun is a
%  kernel function of the form val = kerfun(k, x, nx, y, ny)
%  returning k(s,t) without the speed |z'(t)| factor.
%  s is the segment whose self-interaction is needed.
%  opts.ord controls Alpert quadrature order; other options ignored.
%  Now vectorized for speed! (Most time is spent on kernel evaluation, good)
%
%  If segment.qpblocha is nonempty, segment is treated as grating periodic with
%   Bloch phase qpblocha, and segment parametrization is assumed to spill over
%   correctly into t<0 and t>1. Fixes are done appropriate for nei>0 QP scheme.
%
% Issues:
%  1) Matlab makes local copy of whole of A, which is memory-movement bad!
%  2) segment.Z, Zp, etc eval calls at general t are slow (600us), so code does
%     such requests in bulk (vectorized amortizes object-oriented overhead).
%
% Also see: LAYERPOT.S, LAYERPOT.D, LAYERPOT.T, TEST/TESTLPQUAD

% 1/30/11. Copyright (C) 2011 Alex Barnett
M = size(A,1); N = size(A,2);
if M~=N, warning('Alpert quadr will fail unless matrix square!'); end
qp = ~isempty(s.qpblocha); if qp, a = s.qpblocha; end  % handle grating segs

[tex,wex,nskip] = quadr.QuadLogExtraPtNodes(o.ord); % get log nodes, wei
tex = [-tex(end:-1:1); tex]; wex = [wex(end:-1:1); wex]; % make symmetric
wex = wex/N;      % scale for curve regular node parameter spacing
Na = numel(tex);   % # Alpert pts (both sides)
if N<=2*nskip, warning('not enough regular nodes for Alpert order!'); end
ninterp = o.ord + 3;   % # local regular Lagrange interpolation pts (O'Neil=+2)

% For Alpert, first need to kill off the band-diagonal of A:
% (following loops are over rows, i=targets, j=source density dofs)
if qp  % gratings: kill skipped j's near diag, but *not* SW and NE corners...
  for i=1:N, A(i, [max(i-nskip+1,1):min(i+nskip-1,N)]) = 0; end
  % Now kill SW,NE corner contrib from neighbors that'll be trapezoid summed...
  t = (-nskip+2:0)/N;        % NE corner: t<0 values beyond end of seg
  y = s.Z(t); ny = s.Zn(t); sp = abs(s.Zp(t)); % loc, nor, sp for these src pts
  for i=1:nskip-1
    K = sp(i:end) .* kerfun(k, s.x(i), s.nx(i), y(i:end), ny(i:end));
    jj = (-nskip+1+i:0) + N;  % for this row i, indices j in NE corner
    A(i,jj) = A(i,jj) - a*K/N;   % cancel bloch phase, & weights are all 1/N
  end
  t = 1 + (1:nskip-1)/N;        % SW corner: t>1 values beyond end of seg
  y = s.Z(t); ny = s.Zn(t); sp = abs(s.Zp(t)); % loc, nor, sp for these src pts
  for i=N-nskip+2:N
    jj = 1:i-N+nskip-1; % j inds in SW, also indices in the little t,y, etc list
    K = sp(jj) .* kerfun(k, s.x(i), s.nx(i), y(jj), ny(jj));
    A(i,jj) = A(i,jj) - (1/a)*K/N;   % cancel bloch phase, & weights are all 1/N
  end
  
else  % s closed: kill wrapped band diagonal (ie including SW and NE corners)...
  for i=1:N, A(i, mod(i+[-nskip+1:nskip-1]-1, N)+1) = 0; end
end

% Do set-up for Alpert quadrature extra nodes and weights...
% These arrays are N-by-Na... note that t params can spill over from [0,1]
x = repmat(s.x, [1 Na]); nx = repmat(s.nx, [1 Na]); % target locations & normals
t = repmat(s.t, [1 Na]) + repmat(tex.'/N, [N 1]); % src curve params
if ~qp, t = mod(t,1); end                         % wrap into [0,1]
% Get all curve info from segment at once = fast! (s.Z calls are 600 us, crazy)
y = s.Z(t); ny = s.Zn(t); sp = abs(s.Zp(t));  % speed funcs at src pts
K = sp .* kerfun(k, x, nx, y, ny); % get all kernel evals at once = fast

for l=1:numel(tex)  % loop over Alpert's nonregular quadr pts
  dt = tex(l)/N;    % parameter offset for current quadr pt
  joffs = floor(tex(l)-ninterp/2+1) + (0:ninterp-1); % nearest interp pts
  % interpolation pts and L = row vec of Lagrange basis funcs eval at dt:
  x = joffs/N; w = utils.baryweights(x); L = utils.baryprojs(x, w, dt);
  ii = kron(ones(1,ninterp), 1:N);   % i-indices of A to write over
  jj = ii + kron(joffs, ones(1,N));  % j-indices not yet wrapped into [1,N]
  if qp, wrapphase = a.^-floor((jj-1)/N); end % bloch phase factors of indices
  jj = mod(jj-1, N) + 1;             % j-indices of A to write over
  iii = sub2ind(size(A), ii,jj);     % 1d array indices ditto
  if qp
    A(iii) = A(iii) + wex(l) * wrapphase .* kron(ones(1,ninterp), K(:,l).') .* kron(L, ones(1,N));
  else
    A(iii) = A(iii) + wex(l) * kron(ones(1,ninterp), K(:,l).') .* kron(L, ones(1,N));
  end
end
% Note the affected elements spill out to dist nskip-1+iord/2 from diag!
% (even though the diag band that was set to zero is only to dist nskip-1)
