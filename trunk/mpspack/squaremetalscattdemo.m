% Demonstration of domain-decomp MPS scatt from Dirichlet (sound-soft) square.
% Timo's code testsqscatt2.m, fixed for change in setbc convention, 8/14/08

clear classes
k = 20; % Wavenumber
r = 2.0; % Radius of outer circle    
M = 100; % Number of quadrature points on segments
N=60; % Number of basis fct. in each subdomain
a=1; % Half-Size of the square

% Define all segments
s = segment.polyseglist(M, [1i*r 1i*a a+1i*a a r], 'g');
s=[s(1:4) segment(2*M, [0 r 0 pi/2])];
s = [s rotate(s, pi/2) rotate(s, pi) rotate(s, 3*pi/2)];
s=[s(1:5) s(6:8) s(10) s(11:13) s(15) s(17:18) s(20)];


sdecomp=s([1 6 10 4 5 9 13 16]); % All segments that are part of the domain decomp.
                                 % and not of the square
sq=s([2 8 7 12 11 15 14 3]);     % All segments that are part of the square

% Define the domains
d(1)=domain(s(1:5),[1 1 1 1 1]);
d(2)=domain(s([9 6 7 8 1]),[1 1 1 1 -1]);
d(3)=domain(s([13 10 11 12 6]),[1 1 1 1 -1]);
d(4)=domain(s([16 4 14 15 10]),[1 -1 1 1 -1]);
e = domain([], [], s([5 16 13 9]), -1); 

% Define boundary conditions
sdecomp.setmatch([1 -1],[1 -1]);        % setmatch is now vectorized
sq.setbc(-1,1,0,[]);                    % sq segs will be Dirichlet, vectorized
%        ^--- note sign change here.

% Add basis functions
opts.type='s';
d(1).addnufbbasis(a+1i*a,2/3,-1i,-1-1i,N,k,opts);
d(2).addnufbbasis(-a+1i,2/3,1,1-1i,N,k,opts);
d(3).addnufbbasis(-a-1i,2/3,1i,1+1i,N,k,opts);
d(4).addnufbbasis(a-1i,2/3,-1,-1+1i,N,k,opts);

for dom=d, dom.addregfbbasis(0,N,k); end
e.addmfsbasis(@(z) exp(i*z), -log(1.5*a),2*N,k);

% Setup the problem class
pr=scattering([d,e],[]);
pr.setoverallwavenumber(k);
pr.setincidentwave(-pi/3);

% Solve and print solution
tic; pr.solvecoeffs; fprintf('\tcoeffs done in %.2g sec\n', toc)
fprintf('\tL2 bdry error norm = %g, coeff norm = %g\n', ...
        pr.bcresidualnorm, norm(pr.co))
o.bb=[-3 3 -3 3];
o.dx=0.05;
pr.showthreefields(o);
