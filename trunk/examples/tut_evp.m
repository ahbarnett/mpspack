% Example codes from MPSpack tutorial, also generates EPS figures for this doc
% Barnett 10/8/12
% EIGENVALUE PROBLEMS
clear all classes; verb = 0;

s = segment.smoothnonsym(80, 0.3, 0.2, 3);  % smooth nonsymmm,  Dirichet case
d = domain(s, 1);        % create an interior domain
s.setbc(-1, 'D');        % Dirichlet BC's applied on inside of segment
p = evp(d);              % sets up eigenvalue problem object
d.addlayerpot(s, 'd');          % DLP basis set appropriate for Dir BC
p.solvespectrum([2 10], 'fd');  % Fredholm det method, no modes
p.weylcountcheck(p.kj(1),p.kj,d.perim,d.area);
o.modes = 1; p.solvespectrum([2 10], 'fd', o);  % want modes too
p.showmodes;
if verb, print -depsc2 ../doc/rfnDlow.eps
end

% check L2 norm roughly right: (default is grf)
[uj gx gy di js] = p.showmodes; u = uj(:,:,1); u = u(~isnan(u));
dx = gx(2)-gx(1); sqrt(sum(dx^2*abs(u).^2));
u2 = uj(:,:,2); u2 = u2(~isnan(u2)); sum(dx^2*u2.*u) % crude orthogonality




% Neumann case:
s.setbc(-1, 'N');
d.clearbases;
d.addlayerpot(s, 's');          % SLP set appropriate for Neu BC
o.modes = 1; p.solvespectrum([1 7.8], 'fd', o);
figure;
p.weylcountcheck(p.kj(1),p.kj,-d.perim,d.area); % flip sign of perim for Neu BCs
p.showmodes;
if verb, print -depsc2 ../doc/rfnNlow.eps
end

% check L2 norm roughly right: (default is basis; should be grid-norm'ed)
[uj gx gy di js] = p.showmodes; u = uj(:,:,1); u = u(~isnan(u));
dx = gx(2)-gx(1); sqrt(sum(dx^2*abs(u).^2))
u2 = uj(:,:,2); u2 = u2(~isnan(u2)); sum(dx^2*u2.*u) % crude orthogonality


% high-lying Dirichlet
clear; s = segment.smoothnonsym(80, 0.3, 0.2, 3);
d = domain(s, 1);        % create an interior domain
s.setbc(-1, 'D');        % Dirichlet BC's applied on inside of segment
p = evp(d);              % sets up eigenvalue problem object
d.addlayerpot(s, 'd');          % DLP basis set appropriate for Dir BC
p.updateN(300);
tic; o.eps = 0.1; o.khat = 'r'; o.fhat = 's'; o.modes = 1;
p.solvespectrum([30 31], 'ntd', o);
p.showmodes;
toc
if verb, print -depsc2 ../doc/rfnDk30.eps
end

% Dirichlet triangle
clear; o.kressq = 5;    % corner-packing parameter for Kress reparametrization
s = segment.polyseglist(60, [1, exp(3i*pi/8), exp(5i*pi/4)], 'pc', o);
d = domain(s, 1);        % create an interior domain
s.setbc(-1, 'D');        % Dirichlet BC's applied on inside of segment
p = evp(d);              % sets up eigenvalue problem object
d.addlayerpot(s, 'd');          % DLP basis set appropriate for Dir BC
o.modes = 1; tic; p.solvespectrum([4 20], 'fd', o); toc
p.weylcountcheck(p.kj(1),p.kj,d.perim,d.area);
tic; p.showmodes; toc
if verb, print -depsc2 ../doc/triDlow.eps
end
