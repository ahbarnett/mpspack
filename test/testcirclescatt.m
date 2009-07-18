% Test scattering on a circle

close all
clear all
clear classes

r=1; % Radius of circles
k=10; % Wavenumber

M=300; % Number of points on circle
N=60;  % Number of MFS basis fct. in circle
Rmfs=0.8; % Radius of fundamental solutions inside circle

x0=0; y0=0;
seg=segment(M,[0 r 0 2*pi],'p');
seg.setbc(1,'D',[]);

d=domain([],[],seg,-1);

% Add basis functions

opts.fast=1; opts.eta=k;
Z=@(w) Rmfs*exp(2i*pi*w);
Zp=@(w) 2i*pi*Z(w);
d.addmfsbasis({Z,Zp},N,opts);

pr=scattering(d,[]);
pr.setoverallwavenumber(k);
pr.setincidentwave(-pi/3);



% Solve and print solution
tic; pr.solvecoeffs; fprintf('\tcoeffs done in %.2g sec\n', toc)
fprintf('\tL2 bdry error norm = %g, coeff norm = %g\n', ...
        pr.bcresidualnorm, norm(pr.co));
    

o.bb=[-1.5 1.5 -1.5 1.5];
o.dx=0.01;
o.sepfigs=1;

pr.showthreefields(o);
