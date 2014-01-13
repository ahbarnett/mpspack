%-------------------------------------------
% example of computing the far field using
% MPSPACK.
%
% Copyright (C) 2014 Stuart C. Hawkins
%-------------------------------------------

clear

%-------------------------------------------
% set main parameters
%-------------------------------------------

% wavenumber
kwave = 10;

% parameter that sets separation of sources from the boundary
tau = 5e-2;

% number of points on the boundary at which BC is matched
m = 70;

% number of sources
n = 100;

% set incident wave angle
inc = 0;

%-------------------------------------------
% setup MFS
%-------------------------------------------

% setup structure with the parameters for the MFS method
opts= struct('eta',kwave,'fast',2,'multiplier',2.1,'tau',tau);

% describe the circle in polar coordinates... f gives the radius as a 
% function of the angle and df gives the derivative of the radius
f=@(x) 1.0*ones(size(x));
df=@(x) 1.0*zeros(size(x));

% create boundary segment from the parametrisation
boundary = segment.radialfunc(m, {f,df});

% set a homogeneous Dirichlet boundary condition on the boundary
% ie sound-soft BC
boundary.setbc(1,'D', []);

% setup a domain outside the boundary
d = domain([], [], boundary, -1); 

% setup a MFS basis on the boundary
d.addmfsbasis(boundary, n, opts); 

% initialise the scattering problem
p = scattering(d, []);

% set wavenumber
p.setoverallwavenumber(kwave);

% set incident wave direction
p.setincidentwave(inc);

% solve scattering problem
p.solvecoeffs;

%-------------------------------------------
% plot scattered field
%-------------------------------------------

% setup structure describing the mesh... 
% dx is mesh spacing and bb specifies bounding box
obj = struct('dx',0.08,'bb',3*[-1 1 -1 1]);

% show the solution (scattered wave)
figure(1)
p.showsolution(obj);

%-------------------------------------------
% plot farfield
%-------------------------------------------

clear obj

% setup structure describing the mesh... 
% dx is theta spacing, I think
obj = struct('dx',0.01);

% set whether to plot the real part or imaginary part or neither
% (if neither then the absolute value is plotted)
obj.real = 0;
obj.imag = 0;

% show the far field solution
figure(2)
p.showfarfield(obj);

obj.dB = 1; % show the solution in dB
figure(3); p.showfarfield(obj);

