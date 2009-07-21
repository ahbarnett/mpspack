% SQUAREMETALSCATTDEMO - Exponentially accurate computation of scattering
% from the unit square.

k = 10; % Wavenumber
r = 2.0; % Radius of outer circle    
M = 200; % Number of quadrature points on segments
N=100; % Number of basis fct. in each subdomain
a=1; % Half-Size of the square
rmfs=0.8*r;

% Define all segments
s = segment.polyseglist(M, [1i*r 1i*a a+1i*a a r], 'g');
s=[s(1:4) segment(2*M, [0 r 0 pi/2])];
s = [s rotate(s, pi/2) rotate(s, pi) rotate(s, 3*pi/2)];
s=[s(1:5) s(6:8) s(10:13) s(15) s(17:18) s(20)];


sdecomp=s([1 6 10 4 5 9 13 16]); % All artificial boundaries
sq=s([2 8 7 12 11 15 14 3]);     % All segments that are part of the square
extlist=s([5 16 13 9]);

% Define the domains
d(1)=domain(s(1:5),[1 1 1 1 1]);
d(2)=domain(s([9 6 7 8 1]),[1 1 1 1 -1]);
d(3)=domain(s([13 10 11 12 6]),[1 1 1 1 -1]);
d(4)=domain(s([16 4 14 15 10]),[1 -1 1 1 -1]);
ext = domain([], [], extlist, -1); 

% Define boundary conditions
sdecomp.setmatch([k -k],[1 -1]); % Matching conditions for
                                 % artificial boundaries
extlist.setmatch('diel', 'tm');


% Add basis functions
nuopts.type='s';
d(1).addnufbbasis(a+1i*a,2/3,-1i,-1-1i,N,nuopts);
d(2).addnufbbasis(-a+1i,2/3,1,1-1i,N,nuopts);
d(3).addnufbbasis(-a-1i,2/3,1i,1+1i,N,nuopts);
d(4).addnufbbasis(a-1i,2/3,-1,-1+1i,N,nuopts);

Z=@(t) rmfs*exp(2i*pi*t); Zp=@(t) 2i*pi*Z(t);
opts=struct('eta',k,'fast',1,'nmultiplier',2);
ext.addmfsbasis({Z, Zp},N,opts);


% Setup the problem class
pr=scattering(ext,d);
pr.setoverallwavenumber(k);
pr.setincidentwave(-pi/4);

% Solve and print solution
tic; pr.solvecoeffs; fprintf('\tcoeffs done in %.2g sec\n', toc)
fprintf('\tL2 bdry error norm = %g, coeff norm = %g\n', ...
        pr.bcresidualnorm, norm(pr.co))
o.bb=[-3 3 -3 3];
o.dx=0.05;
pr.showthreefields(o);
