function [x w] = traprule(N)
% TRAPRULE  quadrature points and weights for composite N+1-point trapezoid
% rule

% Copyright (C) 2008, 2009, Alex Barnett, Timo Betcke

x = 2*(0:N)'/N - 1;             
w = (2/N)*ones(1, N+1);
w(1) = 1/N; w(N+1) = w(1);
