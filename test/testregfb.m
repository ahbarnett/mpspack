% Timo's test case which gives NaN on his 64-bit, Matlab 2009a.
% It works fine on my 32-bit Matlab 2008a.

M=16;
N=15;
z0=1+1i;

s=segment(M,[-1.500000000000000e+00 + 7.500000000000001e+00i  -1.949747468305832e+00 + 7.949747468305833e+00i]);


d = domain();
d.addregfbbasis(-1.500000000000000e+00 + 7.500000000000000e+00i, N);
d.k = 10;
d.bas{1}.besselcode = 'r';      % recursivebessel
[A An] = d.bas{1}.eval(s);    % includes endpoint
find(isnan(A))

d.bas{1}.besselcode = 'm';      % validate against Matlab.
[Am Anm] = d.bas{1}.eval(s);
norm(A-Am)
norm(An-Anm)
