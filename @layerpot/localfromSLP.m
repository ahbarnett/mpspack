function A = localfromSLP(s, Jexp)
% LOCALFROMSLP - return matrix taking SLP density values on seg to J (local) exp
%
%  A = localfromSLP(s, Jexp) with s a SLP src segment, and a regfbbasis
%   object Jexp. Jexp contains k, rescale_rad, origin, etc.
%
%  Returns A_{mj} = w_j . sc_m . (i/4) exp(-im.ang(x_j)) H^(1)_m (k|x_j|)
%  for m=-M...M, j=indices of source pts on layer density on segment.
%  sc_m is the J-rescaling factor used in the J-expansion.
%  Note x_j is taken relative to the origin of J-exp
%
%  Requires no layerpot object, but in contrast needs a regFBbasis object.

% Copyright (C) 2008, 2009, Alex Barnett and Timo Betcke


  M = Jexp.N;                       % get max order from regFB object
  sc = Jexp.Jrescalefactors(0:M);   % rescaling factors from regFB obj
  if ~Jexp.real
    sc = [sc(end:-1:2) sc];         % ordering for complex exp basis, -M<=m<=M
  else
    warning 'localfromSLP not implemented for real regFB basis yet!'
    sc = [sc sc(2:end)];            % ordering for real 0:M cos then 1:M sin
  end
  N = numel(s.x);                   % # src pts
  k = Jexp.k;
  x = s.x - Jexp.origin;            % take all coords relative to J-exp origin
  
  % Now evaluate matrix of H^(1)_m (k|x|) . exp(-im.ang(x)) ...
  rotexp = exp(-1i*(1:M)' * angle(x)'); % *=outer-prod, exp(-im.ang(x)), 1..M
  A = [conj(rotexp(end:-1:1,:)); ones(1, N); rotexp];  % stack to m=-M..M
  bes = besselh(repmat([0:M]', [1 N]), repmat(k*abs(x)', [M+1 1]));
  sgn = (-1).^(M:-1:1)';            % alternating signs of negative orders
  A = A .* [bes(end:-1:2,:).*repmat(sgn, [1 N]); bes];  %stack to m=-M..M
  A = A .* repmat((1i/4) * s.w, [2*M+1 1]); % i/4, src quadr wei (no speed!)
  A = A .* repmat(sc.', [1 N]);                 % do J-rescaling
