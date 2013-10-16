function [xm ym info] = gridminfit(f, g, o)
% GRIDMINFIT - use grid to find all minima of min-sing-val-like 1D vector func
%
% This is a global 1d minimizer optimized for vector functions such as the
%   list of  singular or tension values vs k in a MPS-style eigenvalue
%   solver. Zero is special, ie, it is expected minima close to zero in value.
%
% [xm ym] = gridminfit(f, xgrid), where f is handle of a func
%   with interface [ys] = f(x) returning a *vector* of m values (m>=1), and
%   xgrid is a grid, finds all minima of the minimum of the
%   entries in the vector lying in the interval, using an initial search on
%   the grid. iterparabolafit is called on triples of grid samples at local grid
%   minima. If a triple of grid samples has 2nd-smallest entries too close
%   the smallest (based upon maxslope below), recursion onto a smaller uniform
%   grid is done. If f is scalar, no recursion is done.
%   The initial grid need not be uniform, but shouldn't be crazily nonuniform
%   Outputs:
%     xm = row vector of n ordinates of minima
%     ym = m-by-n array of function col-vectors at each of the xm
%
% [xm ym] = gridminfit(f, xgrid, opts) controls various options such as
%      opts.maxslope - function handle to maxslope(x) which returns the
%         max allowable const c in c(x-b) (default 0) 
%      opts.xtol    - desired tolerance in x (default 1e-8). Note this also
%         sets the expected worst-case minima vals as xtol.maxslope.
%      opts.maxit   - max # recursions of forming new grid (default 30)
%      opts.verb    - 0 (silent; default) or 1 (debugging output)
%
% [xm ym info] = ... also returns info struct with info.fevals # func evals
%
% Notes: This supercedes fparamin. Discussion with Taylor Sipple has helped.
% v.1.0 Aug 2011
% v.1.1 added crude degeneracy handling, 10/28/11
% * Could build in a quick check for exact multiplicities based on higher sing
%   vals being less that xtol.maxslope ? Saves lots of recursion
%
% See also: EVP.ITERPARABOLAFIT, EVP.SOLVESPECTRUM

% Copyright (C) 2011, Alex Barnett

if nargin<3 | isempty(o), o = []; end    % process options, defaults...
if ~isfield(o, 'xtol'), o.xtol = 1e-8; end  % default returned arg tolerance
if ~isfield(o, 'maxit'), o.maxit = 30; end
if ~isfield(o, 'verb'), o.verb = 0; end
if ~isfield(o, 'maxslope'), o.maxslope = @(x) 0; end % default not to recurse
if isnumeric(o.maxslope), o.maxslope = @(x) o.maxslope; end % ensure a func!

n = numel(g);
if ~isempty(find(diff(g)<0)), error('grid not everywhere increasing!'); end
if o.verb, fprintf('eval f function at %d grid points...\n', n); end
y = f(g(1)); m = numel(y);     % get m
if o.verb, fprintf('f(%.15g)=%.15g  (min of m=%d entries)\n',g(1),min(y),m); end
yg = nan(m,n); yg(:,1) = sort(y(:));  % preallocate and put in 1st col
% intensive part: eval sorted col vecs from f on grid...
for i=2:n, y = f(g(i)); if o.verb, fprintf('f(%.15g)=%.15g\n',g(i),min(y)); end
  yg(:,i) = sort(y(:)); end    % stack into array
% note, if RAM issue, could keep only lowest 2 values of f at each x ?
info = []; info.ys = yg(1:min(2,m),:); info.xs = g; % keep lowest two sing vals

if o.verb>2, figure; plot(g, yg(1,:), '+-'); % plot lowest
  if m>1, hold on; plot(g, yg(2,:), 'g+'); end, end % plot 2nd lowest

fe = n;   % count fevals
t = zeros(size(g)); % grid pt flags: 1 if used as ctr pt for search, 0 if no
xm = []; ym = []; % locate discrete minima on grid and zoom in or para fit...
for i=1:n
  if i==n, htyp = g(n)-g(n-1); else htyp = g(i+1)-g(i); end % loc grid spacing
  ms = o.maxslope(g(i));
  fm = 1.1*sqrt(htyp^2+o.xtol^2) * ms;   % largest acceptable min at this grid
  if islocmin(i,yg,fm,n)
    if i==1, ii=1:3; elseif i==n, ii=n-2:n;
    else, ii = i-1:i+1; end                          % triple of indices
    if sum(t(ii))==0     % only consider if triple not already flagged
      t(ii(2)) = 1;  % flag the central one
      mingap = nan;  % smallest gap up to 2nd-smallest value (nan is dummy)
      if m>1, mingap = min(yg(2,ii)-yg(1,ii)); end
      if o.verb, fprintf('considering %d %d %d:\n',ii(1),ii(2),ii(3));
        fprintf('htyp=%g, fm=%g, maxslp=%g, mingap=%g\n',htyp,fm,ms,mingap); end
      if o.maxit>0 && m>1 && (g(ii(3))-g(ii(1))>o.xtol) && mingap < 1.1*htyp*ms   % 2nd smallest too close?
        p = o; p.maxit = o.maxit - 1;          % use up one of the recursions
        if o.verb, fprintf('recursing, maxit=%d...\n',o.maxit); end
        [ii t] = fattenlist(ii,t,g,yg,fm); % build out the ii list
        if o.verb, fprintf('\t ii fattened to size %d\n',numel(ii)); end
        shrink = 3;          % factor to shrink the grid by, then recurse...
        [x y in] = evp.gridminfit(f, linspace(g(ii(1)),g(ii(end)),(length(ii)-1)*shrink+1), p);
        fe = fe + in.fevals; info.ys=[info.ys in.ys]; info.xs=[info.xs in.xs];
      else
        p = o; p.maxit = 10;        % for para fit on square of function
        if o.verb, fprintf('%.15g ', [g(ii), yg(1,ii)]); fprintf('\n'); end
        [x y in] = evp.iterparabolafit(@(x) f(x).^2, g(ii), yg(:,ii), p); % get one acc min
        %x,y,in.xs,in.ys,in
        y = sqrt(y); fe = fe + in.fevals; info.xs=[info.xs in.xs];
        info.ys=[info.ys sqrt(in.ys)];% note: sqrt(ys) is kept since sees f^2
        if isempty(y), ndeg=0;
        else
          ndeg = numel(find(y-y(1) < o.xtol*ms)); % estimate degeneracy, crude
        end
        if ndeg>1, x = repmat(x,[1 ndeg]); y = repmat(y,[1 ndeg]); end % dupli!
      end
      if o.verb, fprintf('found %d min(s): ', numel(x));  % plot all pairs...
        if ~isempty(x), a = [x;y(1,:)];
          fprintf('f(%.15g)=%g ',a(:)); fprintf('\n'); end, end
      xm = [xm x]; ym = [ym y];   % stack on end
end, end, end

info.fevals = fe;
return

function a = islocmin(i,yg,fm,n)
if i==1
  a = yg(1,1)<fm && yg(1,1)<yg(1,2);
elseif i==n
  a = yg(1,n)<fm && yg(1,n)<yg(1,n-1);
else
  a = yg(1,i)<yg(1,i-1) && yg(1,i)<yg(1,i+1); % min nr end or i in loc min
end

function [ii t] = fattenlist(ii,t,g,yg,fm) % build out the ii list
% start from ii length 3, grow while values are small enough or bump into ends
% or bump into one next to a t-flagged one.
n = numel(g);
while ii(end)<n && ~t(ii(end)+1) && yg(1,ii(end)+1)<fm % if can grow
  ni = ii(end)+1;   % the new item
  t(ii(end))=1; ii = [ii ni]; % grow to right
  if ni==n-1, ii = [ii n]; t(n-1) = 1; end % fill to right grid end to prevent
  % a 1-grid spacing gap.
end
while ii(1)>1 && ~t(ii(1)-1) && yg(1,ii(1)-1)<fm % if can grow
  ni = ii(1)-1;   % the new item
  t(ii(1))=1; ii = [ni ii]; % grow to left
  if ni==2, ii = [1 ii]; t(2) = 1; end % fill to left grid end
end
