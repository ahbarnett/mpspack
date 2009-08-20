function [x w] = peritrap(N)
% PERITRAP - trapezoid quadrature rule for periodic functions

% Copyright (C) 2008, 2009, Alex Barnett, Timo Betcke

x = 2*(1:N)'/N - 1;
w = 2*ones(1,N)/N;
