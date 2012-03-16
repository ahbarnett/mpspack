% Demo finding (not high) Dir eigenmodes of triangles via Fredholm-det method
% Barnett 1/16/12

clear; N=100; % good up to k=20 @ 1e-12 error w/ q=4 (matrix order is 3N)
o.kressq = 4; % is best for acute triangle Dirichlet interiors
s = segment.polyseglist(N, [1, exp(3i*pi/8), exp(5i*pi/4)], 'pc', o); % give tri
tri = domain(s, 1);
s.setbc(-1, 'D');               % Dirichlet BC's applied on inside: note -1
p = evp(tri);                   % sets up eigenvalue problem object
tri.addlayerpot(s, 'd');        % layerpot basis needed for 'fd' method
o.modes = 1;                    % compute mode representations too
% small-scale calc...
kint = [4 20]; tic; p.solvespectrum(kint, 'fd', o); toc %lowest 25 modes: 46 sec
tic; p.showmodes; toc           % plot the resulting modes

% medium-scale calculation... (around 80 mins)
p.updateN(200);  % good up to k=100 @ 1e-12 error w/ q=4
kint = [4 100]; tic; p.solvespectrum(kint, 'fd', o); toc % errors seem 1e-14 typ
save trik100_modes
tic; o.inds = 762; o.dx = 0.01; p.showmodes(o); toc  % plot highest mode only

if 0 % NtD scaling does not yet work due to needing D, S ops btw multiple segs:
  p.updateN(140); % good up to k=50 @ 1e-12 error w/ q=4
  o = []; o.modes=1; o.eps = 0.1; o.khat='r'; o.fhat='s';
  tic; p.solvespectrum([40 50], 'ntd', o); toc
end

% Neumann case:  only gives 3 digits (o.tol=1e-3) right now
clear; o.kressq = 5; N=130;
s = segment.polyseglist(N, [1, exp(3i*pi/8), exp(5i*pi/4)], 'pc', o); % give tri
tri = domain(s, 1);
s.setbc(-1, 'N');               % Neumann BC's applied on inside: note -1
p = evp(tri);                   % sets up eigenvalue problem object
tri.addlayerpot(s, 's');        % layerpot basis needed for 'fd' method
o.modes = 1;                    % compute mode representations too
% small-scale calc...
kint = [3.5 15]; o.tol=1e-3; tic; p.solvespectrum(kint, 'fd', o); toc
tic; p.showmodes; toc           % plot the resulting modes
