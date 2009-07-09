% PERITRAP - trapezoid quadrature rule for periodic functions

function [x w] = peritrap(N)

x = 2*(1:N)'/N - 1;
w = 2*ones(1,N)/N;
