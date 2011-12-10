function [xm fm info] = iterparabolafit(f, x, y, o)
% ITERPARABOLAFIT - bnded minimize func of 1 var by iterative parabola fit
%
% [xm fm] = iterparabolafit(f, x, y) where f is a function handle to a
%   scalar function with interface y = f(x), and x(1:3) are increasing x-values
%   and y(1:3) corresponding y-values y=f(x), iterates on parabolic fitting
%   to return xm the location of the minimum and fm its value, lying in the
%   interval [x(1) x(3)]. If the minimum appears to lie outside the interval
%   (ie too many rescue attempts), empty sets will be returned.
%
%   If f returns a vector of m values, and y is a m-by-3, the minimum will
%   be used of the entries returned by each call to f,
%   but fm will return the full (sorted) vector of values at the min.
%
% [xm fm] = iterparabolafit(f, x, y, opts) controls various options such as
%      opts.xtol    - desired tolerance in x (default 1e-8)
%      opts.maxit   - max # iterations (default 10)
%      opts.maxresc - max # rescues from outside the interval (default 2)
%      opts.verb    - 0 (silent; default) or 1 (debugging output)
%
% [xm fm info] = ... also returns info struct with info.fevals # func evals,
%      info.ys y values found at info.xs x values computed at.
%
%  Issues:
%   * how handle roundoff in subtraction as rescue ends gets too close?
%    ie add test if b is garbage
%   * allow varargin to pass params to f?
%   * pass out curvature c ?

% Copyright (C) 2011, Alex Barnett

if nargin<4 | isempty(o), o = []; end    % process options, defaults
if ~isfield(o, 'xtol'), o.xtol = 1e-8; end
if ~isfield(o, 'maxit'), o.maxit = 10; end
if ~isfield(o, 'maxresc'), o.maxresc = 2; end
if ~isfield(o, 'verb'), o.verb = 0; end
resc = 0.1;  % factor by which rescued ordinates pulled in from edge

if x(2)<=x(1) || x(3)<=x(2), error('x not increasing!'); end
xoff = x(1); x = x - xoff;   % work relative to left-side point, avoid roundoff
xlo = x(1); xhi = x(3);      % keep original interval
if numel(y)>3, y = sort(y,1); end % sort along the vector axis

info = []; info.xs = []; info.ys = [];
outs = 0; nresc = 0; fe = 0;    % fell-outside flag, # rescues, # fevals
bold = nan;        % keep previous para-fit min estimate
for i=1:o.maxit
  if o.verb, fprintf('x=%.15g %.15g %.15g\n',x(1),x(2),x(3));
    fprintf('y_1=%.15g %.15g %.15g\n',y(1,1),y(1,2),y(1,3));end 
  [a,b,c] = evp.para_fit(x, y(1,:));
  if o.verb, fprintf('para_fit a,b,c = %.15g  %.15g %.15g\n',a,b,c); end
  if abs(b-bold)<o.xtol, yb = sort(f(b + xoff)); yb = yb(:); % col vec
    fe = fe + 1; info.xs=[info.xs,b+xoff]; info.ys=[info.ys yb(1:2)];
    break; end % two close b's, done
  if b<xlo, if o.verb, fprintf('b fell off left\n', b); end
    if nresc<o.maxresc
      b = (1-resc)*xlo + resc*x(2); nresc = nresc+1; % rescue back insid
    else, outs = 1; break; end
  elseif b>xhi, if o.verb, fprintf('b fell off right\n', b); end
  if nresc<o.maxresc
    b = (1-resc)*xhi + resc*x(2); nresc = nresc+1; % rescue
  else, outs = 1; break; end
  end
  if min(abs(b-x))<o.xtol/2, b = b + o.xtol/2; end  % jiggle if too close
  bold = b;
  yb = sort(f(b + xoff)); yb = yb(:); % col vec
  fe = fe + 1; info.xs=[info.xs,b+xoff]; info.ys=[info.ys yb(1:2)];
  if o.verb, fprintf('f_1(%.15g)=%.15g\n',b,yb(1)); end
  [dummy j] = sort(abs(x-b)); j = j(1:2); % indices of two old x's closest to b
  x = [x(j) b]; y = [y(:,j) yb];  % add the new point
  [x j] = sort(x); y = y(:,j);    % resort to ascending ordinates
end

if ~outs
  x(4) = b; y(:,4) = yb;
  [dummy j] = sort(y(1,:)); % sort on just lowest entries!
  y = y(:,j); x = x(j); xm = x(1); fm = y(:,1);  % keep the best point
  %xm = b; fm = yb;  % keep last fit
else     % minima not in interval
  xm = []; fm = [];
end

xm = xm + xoff;  % restore the offset

info.fevals = fe;
