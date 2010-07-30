function h = arrow(x, y, varargin)
% ARROW - draw a 2D arrow from (x1,y1) to (x2,y2).
%
% h = arrow(x, y, ...) where x and y are 2-element arrays puts head at the
%  second element.

z = x+1i*y;
d = x(2)-x(1) + 1i*(y(2)-y(1));
l = 0.2;                        % size of head
h = plot([z(1) z(2)+l*d*[0 -1+.7i, 0 , -1-.7i]], varargin{:});
