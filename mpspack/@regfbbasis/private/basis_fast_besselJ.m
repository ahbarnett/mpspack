function [phi phin] = basis_fast_besselJ(k, b, p, pn)
% function [phi phin] = basis_fast_besselJ(k, b, p, pn)
%
% Simultaneous eval of J-bessel basis (local expansion) at points p,
% and its pn-directional derivative. b gives basis struct, k wavenumber.
% Outputs: phi, phin are size M-by-N (M=#pts, N=#basis funcs).
% Test this code against slower basis.m with test_basis_fast_besselJ.m
% Now uses fast exp(1i...) via repeated rotations, so bessel&exp comparable
%
% barnett 2/29/08, fixed r=0 case causing NaN 4/3/08

if size(p,1)~=2                         % sanity checks
  error('p must be 2-by-something!');
end
if nargin>3
  if size(pn)~=size(p)
    error('pn must be same size as p!');
  end
end
M = size(p,2);                          % how many locations to evaluate at

N = b.N;
n = -b.maxorder + (0:N-1);              % set of desired orders as row vec
MJ = max(abs(n)) + 1;                   % largest order including +1 for Jderiv
sc = ones(1,N);
if isfield(b, 'rescale') & b.rescale    % scale factors for J basis funcs
  sc = 1 ./ Jscalefactor(n, k); 
end
%max(sc)
r = sqrt(sum(p.^2));                    % set of distances
J = fast_besselJ_recurrence(k*r, MJ);   % get all bessels fast
nm = repmat(n, [M 1]);                  % const-cols matrix of n's
if 0       % conventional naive exp
  ang = exp(1i*nm .* repmat(angle(p(1,:) + 1i*p(2,:)).', ...
                            [1 N]));    % ang func (eaty) at all n
else       % intelligent exps using fact that everything is complex rotations
  ang = fast_angle_exp(n, p(1,:) + 1i*p(2,:));
end
idx = MJ-b.maxorder+(1:N);              % indices of n's within list -MJ:MJ
phi = repmat(sc, [M 1]) .* ang .* J(:,idx);    % basis values
if nargout>1                            % eval directional deriv
  i = find(r==0); r(i) = 1;             % zero-r indices, use dummy r val
  invr = 1./r;                          % since dividing slow
  invrm = repmat(invr.', [1 N]);        % const-rows matrix of same
  nr = sum(pn.*p).*invr;
  nt = sum(pn.*[p(2,:); -p(1,:)]).*invr;  % (r,t) components of pn
  phin = repmat(sc, [M 1]) .* ang .* (...
      repmat(nr.', [1 N]) .* (k*J(:,idx-1) - nm.*J(:,idx).*invrm)...
      - repmat(nt.', [1 N]) .* invrm .* J(:,idx) .* 1i.*nm);
  % note above how n=0 not handled specially... unlike in basis.m
  if ~isempty(i)                        % explicitly overwrite the r=0 cases
    phin(i,:) = 0;                      % |n|=1 use slope of J_1, other n's 0
    if n(end)>0, phin(i,b.maxorder+2) = k/2*(pn(1,i)+1i*pn(2,i)).'; end % n = +1
    if n(1)<0, phin(i,b.maxorder) = -k/2*(pn(1,i)-1i*pn(2,i)).'; end    % n = -1
  end  % the above formulae for r=0 found by trial & error via test_evalbasis.m
end
           

