function Rjn = kress_Rjn(n)
% function Rjn = kress_Rjn(n)
%
% return 2n length vector of R_j^(n) for j=0...2n-1. Takes O(n ln n) work
% and O(n) storage, using fft for the trig sum.
% Note the full R_{|i-j|}^(n) matrix is then circulant(kress_Rjn(N/2)).
% See Kress MCM 1991 paper or Lin Int Eqn book p. 210
% barnett 2/6/08

if mod(2*n,2)==1, disp('kress_Rjn: N=2n must be even!'); return; end

m = 1:n-1;
Rjn = -2*pi*ifft([0 1./m 1/n 1./m(end:-1:1)]);

% old direct O(n^2) code
%j = 0:2*n-1;
%Rjn = cos(j*pi)/n;
%for m=1:n-1
%  Rjn = Rjn + 2*cos(m*j*pi/n)/m;
%end
%Rjn = -pi*Rjn/n;
