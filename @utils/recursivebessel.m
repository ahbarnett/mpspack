function J = recursivebessel(M,x)
% Computes J bessel at all orders 0<=M for multiple arguments x, using
% Miller's method of downwards-stable recurrence relations (Num Rec sec 5.4,
% and 6.5, 3rd Ed). Recurrence is started so absolute errors are around 1e-16.
% (Relative errors will be unbounded in evanescent region).
% Limited to |x|, M < 1e8 or some huge number. Trick to avoid underflow is by
% for each x starting the true recurrence only when needed, but before this
% a fake O(1) recurrence is used.
%
% output: J is N-by(M+1), where N is numel(x)     [note changed from Alex code]
%
% barnett 2/28/08
% Modified for inclusion in mpspack by Timo Betcke 13/07/08
% debugged nst check at end, barnett 7/27/08

N = numel(x);
x = reshape(x, [N 1]);                          % make col vec
d = max(abs(x));
eps = 1e-16;               % roughly e_mach, beneath which J_n=1 (n=0) or 1.
nz = find(abs(x) > eps);         % indices of effectively nonzero args
z = find(abs(x) < eps);
x(z) = 1;                  % for these dummy values J is overwritten at end
nrs = 10;          % how often to rescale, allows n/x<1e20, w/ 1e308 overflow
% Airy asym exp in trans regions, A&S 9.3.23, Ai(2^(1/3)*z) < 1e-16 for z=12...
% ...but I bumped it up to z=18 since need rel acc in evan region (rescaled!)
ncutoff = @(x) ceil(x + 18 * (x+1).^(1/3));  % irrelevant for orders beyond this
nst = ncutoff(d);                               % starting order (using Airy)
ncutoffx = ncutoff(x);                          % starting orders for each x
noff = 1;                                       % order index n offset
% in the following, unity is used to initialize all recurrence values...
J = ones(N, noff+nst);                          % store orders 0,1,..,nst
invx = 1./x;                                    % only do the divide once
J(:,noff+nst-1) = 2*nst*J(:,noff+nst).*invx;    % start using J_{nst+1} = 0
for no=nst-1:-nrs:1                             % downwards J stable recurrence
  nl = min(nrs, no);                            % how many to do in inner loop
  % Note the hack to add 5 to ncutoffx gives extra acc but don't make eg 10
  % since overflows in the evan region at small x...
  i = find(no-nl-5 > ncutoffx);                 % indices of x vals not to recur
  invxreg = invx; invxreg(i) = 1./(no-nl/2);    % n-`regularized' (safe) 1/x
  for m=0:nl-1                                  % recur nrs times w/ safe update
    n = no - m;
    J(:,noff+n-1) = 2*n*J(:,noff+n).*invxreg - J(:,noff+n+1);
    % Note that J(:,...) was much faster than indexed filling J(i,...)
  end
end
norms = J(:,noff) + 2*sum(J(:,noff+(2:2:nst)), 2); % compute known unity sum
% normalize using fact that n-sums should give 1...
J = J .* repmat(1./norms, [1 size(J,2)]);
% reinsert the effectively zero args...
J(z,:) = repmat([1 zeros(1,nst)], [numel(z) 1]);

% now trim depending on how many orders needed... (Alex reinserted nst check)
if M<=nst
   J = J(:,1:noff+M);
else
  extran = M-nst;                                % how many cols to add in n
  J = [J(:,1:noff+nst) zeros(N, extran)];
end
