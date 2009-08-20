% Test regular Fourier-Bessel basis. Timo's early code.
% Fixed up, Barnett 7/16/09. Also see cases 1:3 in testbasis.m

k=2;
N=5;

% Create a bunch of random points plus normal direction

M=10;
z=randn(M,1)+1i*randn(M,1);
nz=randn(M,1)+1i*randn(M,1); nz=nz./abs(nz);
pts=pointset(z,nz);

% Create the matrix of Fourier-Bessel fct. up to order 
R=abs(z); T=angle(z);
A=[besselj(0,k*R),besselj(1:N,k*R).*cos(T*(1:N)),besselj(1:N,k*R).*sin(T*(1:N))];

% Now perturb z slightly in the x and y directions to get approximate
% derivatives

zx=z+1E-8*ones(size(z));
zy=z+1i*1E-8*ones(size(z));
Rx=abs(zx); Tx=angle(zx);
Ax=[besselj(0,k*Rx),besselj(1:N,k*Rx).*cos(Tx*(1:N)),besselj(1:N,k*Rx).*sin(Tx*(1:N))];
Ax=(Ax-A)/(1E-8);
Ry=abs(zy); Ty=angle(zy);
Ay=[besselj(0,k*Ry),besselj(1:N,k*Ry).*cos(Ty*(1:N)),besselj(1:N,k*Ry).*sin(Ty*(1:N))];
Ay=(Ay-A)/(1E-8);
nx1=real(nz); nx2=imag(nz);
An=repmat(nx1,1,2*N+1).*Ax+repmat(nx2,1,2*N+1).*Ay;

% Compare with output from regfbbasis.eval
% The following edited by Alex for updated interface, 7/16/09:
opts.real=1;
opts.besselcode='r';
b = regfbbasis(0,N,opts); b.doms = domain(); b.doms.k = k;
[AA AAx AAy]=b.eval(pts);
% I assume this is what you wanted...
norm(A-AA), norm(Ax-AAx), norm(Ay-AAy)
