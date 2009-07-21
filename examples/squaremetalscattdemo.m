% SQUAREMETALSCATTDEMO - Exponentially accurate computation of
% time-harmonic scattering on the unit square.

% Parameters of the problem
% -------------------------

k = 50;      % Wavenumber
r = 1.0;     % Radius of outer circle    
M = 200;     % Number of quadrature points on segments
N=100;       % Number of basis fct. in each subdomain
a=.5;        % Half-Size of the square
rmfs=0.8*r;  % Radius of the fundamental solutions curve

% Definition of the mesh
% ----------------------

% Define all segments

s = segment.polyseglist(M, [1i*r 1i*a a+1i*a a r], 'g');
s=[s(1:3) segment(2*M, [0 r 0 pi/2])];
s = [s rotate(s, pi/2) rotate(s, pi) rotate(s, 3*pi/2)];


sdecomp=s([1 4 5 8 9 12 13 16]); % All artificial boundaries
extlist=s([4 8 12 16]);          % Segments forming the outer circle

% Define the domains
d=domain.empty(4,0);
for j=1:4, d(j)=domain(s(1+mod(4*(j-1)+[0 1 2 12 3],16)),[1 1 1 -1 1]); end
ext = domain([], [],extlist(end:-1:1), -1); 
    
% Boundary conditions between elements
% ------------------------------------

sdecomp.setmatch([k -k],[1 -1]);


% Define the basis functions
% --------------------------

% Fractional Bessel functions
nuopts=struct('type','s','rescale_rad',1);
for j=1:4, 
    d(j).addnufbbasis((1i)^(j-1)*(a+1i*a),2/3,(1i)^(j-1)*(-1i),(1i)^(j-1)*(-1-1i),N,nuopts);
end

% 
% d(1).addnufbbasis(a+1i*a,2/3,-1i,-1-1i,N,nuopts);
% d(2).addnufbbasis(-a+1i*a,2/3,1,1-1i,N,nuopts);
% d(3).addnufbbasis(-a-1i*a,2/3,1i,1+1i,N,nuopts);
% d(4).addnufbbasis(a-1i*a,2/3,-1,-1+1i,N,nuopts);

% Fundamental Solutions
Z=@(t) rmfs*exp(2i*pi*t); Zp=@(t) 2i*pi*rmfs*exp(2i*pi*t);
opts=struct('eta','k','fast',1,'nmultiplier',2);
ext.addmfsbasis({Z, Zp},N,opts);


% Setup the problem class
% -----------------------

pr=scattering(ext,d);
pr.setoverallwavenumber(k);
pr.setincidentwave(-pi/4);

% Solve and plot solution
% -----------------------

tic; pr.solvecoeffs; fprintf('\tcoeffs done in %.2g sec\n', toc)
fprintf('\tL2 bdry error norm = %g, coeff norm = %g\n', ...
        pr.bcresidualnorm, norm(pr.co))
o.bb=[-1.5 1.5 -1.5 1.5];
o.dx=0.01;

[ui gx gy] = pr.gridincidentwave(o);
us = pr.gridsolution(o);

figure;
imagesc(gx, gy, real(ui+us)); title('Full Field (Real Part)');
c = caxis; caxis([-1 1]*max(c));
axis equal tight;
colorbar;
set(gca,'ydir','normal'); 
