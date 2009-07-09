% TRAPRULE  quadrature points and weights for composite N+1-point trapezoid rule

function [x w] = traprule(N)

x = 2*(0:N)'/N - 1;             
w = (2/N)*ones(1, N+1);
w(1) = 1/N; w(N+1) = w(1);
