% Test scattering on a circle

close all

r=1; % Radius of circles
k=10; % Wavenumber

M=300; % Number of points on circle
NN=1:60;  % List for the numbers of basis functions
Rmfs=0.8; % Radius of fundamental solutions inside circle

x0=0; y0=0;
seg=segment(M,[0 r 0 2*pi],'p');
seg.setbc(1,'D',[]);

d=domain([],[],seg,-1);

% Add basis functions

N=1;

opts.fast=1; opts.eta=k;
Z=@(w) Rmfs*exp(2i*pi*w);
Zp=@(w) 2i*pi*Z(w);
d.addmfsbasis({Z,Zp},N,opts);

pr=scattering(d,[]);
pr.setoverallwavenumber(k);
pr.setincidentwave(-pi/3);


resvec=zeros(size(NN)); % Vector with residuals

for j=1:length(NN),
    % Loop that goes through different numbers for N and computes the
    % corresponding residual
    
    % Update the number of basis functions
    pr.updateN(NN(j));
    
    % Solve and print solution
    tic; pr.solvecoeffs; fprintf('\tcoeffs done in %.2g sec for N=%i\n', toc,NN(j))
    resvec(j)=pr.bcresidualnorm;
    
    fprintf('\tL2 bdry error norm = %g, coeff norm = %g\n', ...
        resvec(j), norm(pr.co));
    
end
    

o.bb=[-1.5 1.5 -1.5 1.5];
o.dx=0.01;
o.sepfigs=0;

pr.showthreefields(o);

figure; semilogy(NN,resvec,'k-'); hold on;
convrate=Rmfs/r;
semilogy(NN,convrate.^(NN),'r--');
legend('Measured convergence','Expected convergence');
