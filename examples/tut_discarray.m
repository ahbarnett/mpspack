% Scattering from 3x3 array of discs, as in Tutorial section 10.
% pasted from tutorial.pdf by Barnett 8/2/17

N1=3; N2=3; % Number of scatterers in each direction
r=1; % Radii of circles
a=3; % Distance of midpoints of neighboring circles in each dimension
k=10; % Wavenumber
M=300; % Number of points on each circle
N=150; % Number of MFS basis fct. in each circle
Rmfs=0.8*r; % Radius of fundamental solutions inside circles

y0=0;
s=segment.empty(N1*N2,0);
for i=1:N1,
  x0=0;
  for j=1:N2
    seg=segment(M,[x0+1i*y0 r 0 2*pi],'p');
    seg.setbc(1,'D',[]);
    s((i-1)*N2+j)=seg;
    x0=x0+a;
  end
  y0=y0+a;
end

seg.setbc(1,'D',[]);
o.normals=0;
plot(s,1,o);

d=domain([],[],num2cell(s),num2cell(-1*ones(N1*N2,1)));

x0=0; y0=0; opts.fast=1;
for i=1:N1,
  x0=0;
  for j=1:N2
    Z=@(w) Rmfs*exp(2i*pi*w)+x0+1i*y0;
    Zp=@(w) Rmfs*2i*pi*exp(2i*pi*w);
    d.addmfsbasis({Z,Zp},N,opts);
    x0=x0+a;
  end
  y0=y0+a;
end

pr=scattering(d,[]);
pr.setoverallwavenumber(k);
pr.setincidentwave(-pi/3);

tic; pr.solvecoeffs; fprintf('\tcoeffs done in %.2g sec\n', toc)
fprintf('\tL2 bdry error norm = %g, coeff norm = %g\n', ...
pr.bcresidualnorm, norm(pr.co));

o.bb=[-a N2*a -a N1*a];
o.dx=0.05;
o.sepfigs=1;
pr.showthreefields(o);
