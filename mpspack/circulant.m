function A = circulant(x)
% function A = circulant(x)
%
% return square circulant matrix with first row x
% barnett 2/5/08

x = x(:);
A = toeplitz([x(1); x(end:-1:2)], x);
