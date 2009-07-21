function A = localfromDLP(s, Jexp)
% LOCALFROMDLP - return matrix taking DLP density values on seg to J (local) exp
%
%  A = localfromDLP(s, Jexp) with s a DLP src segment, and a regfbbasis
%   object Jexp. Jexp contains k, rescale_rad, origin, etc.
%
%  Returns A_{mj} = w_j . sc_m . (ik/8) .
%                  [exp(i(1-m).ang(x_j)) H^(1)_{m-1} (k|x_j|) exp(-i.ang(n_j))
%                 - exp(i(-1-m).ang(x_j)) H^(1)_{m+1} (k|x_j|) exp(i.ang(n_j))]
%  for m=-M...M, j=indices of source pts on layer density on segment.
%  sc_m is the J-rescaling factor used in the J-expansion.
%  Note x_j is taken relative to the origin of J-exp
%
%  Requires no layerpot object, but in contrast needs a regFBbasis object.

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
  
  % Now evaluate rotexp and bes matrices which have orders m=-M-1...M+1 :
  rotexp = exp(-1i*(1:M+1)' * angle(x)');% *=outer-prod, exp(-im.ang(x)), 1..M
  rotexp = [conj(rotexp(end:-1:1,:)); ones(1, N); rotexp]; % now m=-M-1..M+1
  bes = besselh(repmat([0:M+1]', [1 N]), repmat(k*abs(x)', [M+2 1]));%0..M+1
  sgn = (-1).^(M+1:-1:1)';            % alternating signs of negative orders
  bes = [bes(end:-1:2,:).*repmat(sgn, [1 N]); bes]; % now m=-M-1..M+1
    
  A = bes(1:end-2,:).*rotexp(1:end-2,:).*repmat(exp(-1i*angle(s.nx).'), [2*M+1 1]) - bes(3:end,:).*rotexp(3:end,:).*repmat(exp(1i*angle(s.nx).'), [2*M+1 1]);
  A = A .* repmat((1i*k/8) * s.w, [2*M+1 1]); % ik/8, src quadr wei (no speed!)
  A = A .* repmat(sc.', [1 N]);                 % do J-rescaling
