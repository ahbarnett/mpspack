function A = localfromSLP(s, Jexp)
% LOCALFROMSLP - return matrix taking SLP density values on seg to local (J) exp
%
%  A = localfromSLP(s, Jexp) with s a SLP src segment, and a regfbbasis
%   object Jexp. Jexp contains k, rescale_rad, origin, etc.
%
%  Returns, for k=0,
%   A_{0j} = w_j . sc_0 . (1/2.pi) log 1/|x_j|
%   A_{mj} = w_j . sc_m . (1/4.pi.m) (-x_j)^{-m},   m>0
%   A_{-m,j} = (-1)^m . conj(A_{mj}),   m>0
%  or otherwise,
%    A_{mj} = w_j . sc_m . (i/4) exp(-im.ang(x_j)) H^(1)_m (k|x_j|)
%  for m=-M...M, j=indices of source pts on layer density on segment.
%  sc_m is the J-rescaling factors used in the J-expansion.
%  Note x_j is taken relative to the origin of J-exp
%
%  Requires no layerpot object, but in contrast needs a regFBbasis object.
%
% See also: test/testlayerpotJfilter.m

% Copyright (C) 2008 - 2012, Alex Barnett and Timo Betcke

% Laplace case added 3/16/12

  M = Jexp.N;                       % get max order from regFB object
  N = numel(s.x);                   % # src pts
  k = Jexp.k;
  x = s.x - Jexp.origin;            % take all coords relative to J-exp origin
  
  if k==0   % Laplace case, like a M2L, the monopole term only, then Re part...
    A = [-2*log(abs(x)') ; (repmat(x.', [M 1]) .^ repmat(-(1:M)', [1 N])) .* repmat(1./(1:M)', [1 N]) ];
    sgn = (-1).^(M:-1:1)';          % signs same as Helmholtz case.
    A = [conj(A(end:-1:2,:)).*repmat(sgn, [1 N]); A];   % stack to m=-M..M
    prefac = 1/(4*pi);               % prefactor (ncludes 1/2 for Re part)
  else   % Now evaluate matrix of H^(1)_m (k|x|) . exp(-im.ang(x)) ...
    rotexp = exp(-1i*(1:M)' * angle(x)'); % *=outer-prod, exp(-im.ang(x)), 1..M
    A = [conj(rotexp(end:-1:1,:)); ones(1, N); rotexp];  % stack to m=-M..M
    bes = besselh(repmat((0:M)', [1 N]), repmat(k*abs(x)', [M+1 1]));
    prefac = 1i/4;
    sgn = (-1).^(M:-1:1)';            % alternating signs of negative orders
    A = A .* [bes(end:-1:2,:).*repmat(sgn, [1 N]); bes];  % stack to m=-M..M
  end
  A = A .* repmat(prefac * s.w, [2*M+1 1]);     % src quadr wei (no speed!)

  if Jexp.rescale_rad~=0
    sc = Jexp.Jrescalefactors(0:M);   % rescaling factors from regFB obj
    if ~Jexp.real
      sc = [sc(end:-1:2) sc];         % ordering for complex exp basis, -M<=m<=M
    else
      warning 'localfromDLP not implemented for real regFB basis yet!'
      sc = [sc sc(2:end)];            % ordering for real 0:M cos then 1:M sin
    end
    A = A .* repmat(sc.', [1 N]);                 % do J-rescaling
  end
