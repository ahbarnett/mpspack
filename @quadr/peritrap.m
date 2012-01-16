function [x w] = peritrap(N)
% PERITRAP - trapezoid quadrature rule for periodic functions

% Copyright (C) 2008 - 2012, Alex Barnett, Timo Betcke

%x = 2*(1:N)'/N - 1;  % original, unsymmetric
x = 2*(1:N)'/N - 1 - 1/N; % shifts half a grid-point back, symmetric about .5
w = 2*ones(1,N)/N;
