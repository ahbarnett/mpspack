  function [x,w] = clencurt(N)
% CLENCURT  nodes x (Chebyshev points) and weights w
%           for Clenshaw-Curtis quadrature. FFT version.
%
%   (was: Trefethen book, modified by Barnett to return x in increasing order)
%
%   Now: FFT version using Fourier series for |sin(theta)|, by Barnett

% Copyright (C) 2008, 2009, 2010, Alex Barnett, Timo Betcke

theta = pi*(N:-1:0)'/N; x = cos(theta); % note order opposite to Trefethen
W = kron(-1./((1:floor(N/2)).^2-1/4), [0 1]);   % works for even or odd
if mod(N,2)==1, W = [W 0]; end  % include extra pi-freq term if odd
w = ifft([4 W W(end-1:-1:1)]);  % 4 is the zero-freq term
w = [w(1)/2 w(2:N) w(1)/2];     % endpoints get 1/2 weight since want 1/2 circle
