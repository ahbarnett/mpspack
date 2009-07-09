
%
% ... construct the trapezoidal quadrature formula and apply end
% correction formula and verify Kapur-Rokhlin's corrections.
%

n=100; m=10;
[x,w]=kapurtrap(n,m);
%

sum(w)-1

%
nmax=m-1; 
rints_exact=1./(1:(nmax+1))';

%% Test non-periodic functions

%
% ... integrate polynomials up to degree nmax
%
rints=zeros(nmax+1,1);
for i=0:nmax
  f=x.^i;
  rints(i+1)=sum(w.*f);
end
errs=rints-rints_exact

%
% ... integrate exp(x)
%
rint_exact=exp(1)-exp(0);
f=exp(x);
rint=sum(w.*f);
err=rint-rint_exact

%% Test periodic functions

a=0;
b=2*pi;

x=x*(b-a)+a;
w=w*(b-a);

rint_exact=0;
f=cos(x);
rint=sum(w.*f);
err=rint-rint_exact
