function D = perispecdiffrow(N)
% PERISPECDIFFROW - row 1 of N-pt spectral 2pi-periodic differentiation matrix

% Copyright (C) 2008, 2009, Alex Barnett, Timo Betcke
if mod(N,2)==1, error 'perispecdiffrow: N must be even!', end
tj = 2*pi/N*(0:N-1);
D = (-1).^(1:N) .* cot(tj/2) / 2;   % note overall - sgn due to 1st row not col
D(1) = 0;              % kill the Inf
