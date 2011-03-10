function h = arrow(x, y, varargin)
% ARROW - draw a 2D arrow from (x1,y1) to (x2,y2).
%
% h = arrow(x, y) where x and y are 2-element arrays puts head at the
%  second element.
%
% h = arrow(x, y, P1, P2, ....) passes params to plot line drawing command.
%
% h = arrow(x, y, opts, P1, P2, ...) passes P1, P2, etc to plot command, but
%  where opts is a struct, also controls the following options:
%  opts.headalong : fraction of the way along to place the head (default 1).
%  opts.headsize : fractional size of head (default 0.2)

% Copyright (C) 2010 - 2011, Alex Barnett
z = x+1i*y;
d = x(2)-x(1) + 1i*(y(2)-y(1));
l = 0.2;                        % size of head
headalong = 1.0;                % location of head, default
argptr = 1; opts  = [];
if numel(varargin)>0 & isstruct(varargin{1})
  argptr = 2; opts = varargin{1}; end
if ~isempty(opts) & isfield(opts,'headsize'), l = opts.headsize; end
if ~isempty(opts) & isfield(opts,'headalong'), headalong = opts.headalong; end
h = plot([z(1) z(2) z(1)+(z(2)-z(1))*headalong+l*d*[0 -1+.7i, 0 , -1-.7i]], varargin{argptr:end});
