function A = localfromDLP(s, Jexp, o)
% LOCALFROMDLP - S2L matrix taking DLP density values on seg to local (J) exp
%
%  A = localfromDLP(s, Jexp) with s a DLP src segment, and a regfbbasis
%   object Jexp. Jexp contains k, rescale_rad, origin, etc.
%
%  Returns, for k=0,
%   ...(to do)...
% 
% or otherwise, A_{mj} = w_j . sc_m . (ik/8) .
%                  [exp(i(1-m).ang(x_j)) H^(1)_{m-1} (k|x_j|) exp(-i.ang(n_j))
%                 - exp(i(-1-m).ang(x_j)) H^(1)_{m+1} (k|x_j|) exp(i.ang(n_j))]
%  for m=-M...M, j=indices of source pts on layer density on segment.
%  sc_m is the J-rescaling factor used in the J-expansion.
%  Note x_j is taken relative to the origin of J-exp
%
%  A = localfromDLP(s, Jexp, o) controls options:
%      o.fast = 0,1,2: chooses Matlab's Hankel, or Rokhlin codes.
%
%  Requires no layerpot object, but in contrast needs a regFBbasis object.
%
% Note: formulae look different due to 
%
% See also: test/testlayerpotJfilter.m

% Copyright (C) 2008 - 2012, Alex Barnett and Timo Betcke

% Laplace case added 3/16/12
% Hankel recurrence (no overflow protect) 11/8/12

  M = Jexp.N;                       % get max order from regFB object
  N = numel(s.x);                   % # src pts
  k = Jexp.k;
  x = s.x - Jexp.origin;            % take all coords relative to J-exp origin
  if nargin<3 | ~isfield(o,'fast'), o.fast = 0; end
  
  if k==0
    A = repmat(x.', [M+1 1]) .^ repmat(-(1:M+1)', [1 N]) .* repmat((-1/(4*pi))*s.w .* s.nx.', [M+1 1]);  % prefactor (ncludes 1/2 for Re part)
    A(1,:) = 2*real(A(1,:));          % m=0 term
    sgn = (-1).^(M:-1:1)';    % signs same as Helmholtz case.
    A = [conj(A(end:-1:2,:)).*repmat(sgn, [1 N]); A];   % stack to m=-M..M
  else  
    % Now evaluate rotexp and bes matrices which have orders m=-M-1...M+1 :
    rotexp = exp(-1i*(1:M+1)' * angle(x)');%*=outer-prod, exp(-im.ang(x)), 1..M
    %rotexp = zeros(M+2,numel(x)); rotexp(1,:) = 1; % recurrence for exps...
    %rotexp(2,:) = exp(-1i * angle(x)');
    %for n=1:M, rotexp(n+2,:) = rotexp(n+1,:).*rotexp(n+1,:)./rotexp(n,:); end
    %R = [conj(rotexp(end:-1:2,:)); rotexp];  % stack to m=-M-1..M+1
    R = [conj(rotexp(end:-1:1,:)); ones(1, N); rotexp]; % now m=-M-1..M+1
    %bes = besselh(repmat([0:M+1]', [1 N]), repmat(k*abs(x)', [M+2 1]));%0..M+1
    bes = zeros(M+2,N);
    if o.fast==2   % now just eval orders 0,1
      [bes(1,:) bes(2,:)] = utils.greengardrokhlinhank106(k*abs(x)');
    elseif o.fast==1
      [bes(1,:) bes(2,:)] = utils.greengardrokhlinhank103(k*abs(x)');
    else
      bes(1:2,:) = besselh(repmat([0;1], [1 N]), repmat(k*abs(x)', [2 1]));
    end
    tinvz = 2./(k*abs(x)');  % row vec of 2/z
    for n=1:M     % recurrence (stable upwards I think!) up to order M+1
      bes(n+2,:) = n*tinvz.*bes(n+1,:) - bes(n,:);
    end
    sgn = (-1).^(M+1:-1:1)';            % alternating signs of negative orders
    bes = [bes(end:-1:2,:).*repmat(sgn, [1 N]); bes]; % now m=-M-1..M+1
    
    A = bes(1:end-2,:).*R(1:end-2,:).*repmat(exp(-1i*angle(s.nx).'), [2*M+1 1]) - bes(3:end,:).*R(3:end,:).*repmat(exp(1i*angle(s.nx).'), [2*M+1 1]);
    A = A .* repmat((1i*k/8) * s.w, [2*M+1 1]); % src quadr wei (no speed!)
  end
  
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
