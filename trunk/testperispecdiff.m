% test periodic spectral differentiation matrix
% Barnett 1/21/09

clear all classes
N = 50;
tj = 2*pi/N*(1:N)';
f = sin(3*tj);    fp = 3*cos(3*tj);         % trial function
D = circulant(quadr.perispecdiffrow(N));
%figure; imagesc(D);
%figure; plot(tj, [fp D*f], '+-');
norm(D*f - fp)

