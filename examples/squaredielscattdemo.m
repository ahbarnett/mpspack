% use scattering class to do Timo's square scattering in 1 page, barnett 7/28/08

clear all
k = 5;
R = 1.0;     % radius of enclosing circle; square will be side 1 centered at 0
M = 25;    % ...... set up segments
s = segment.polyseglist(M, [.5 .5+.5i .5i 1i*R], 'g');
a = segment(2*M, [0 R pi/2 0]);              % CW arc from 12 to 3 o'clock
s = [s(1:3) a];                              % repeating unit
s = [s rotate(s, pi/2) rotate(s, pi) rotate(s, 3*pi/2)]; % NB: not pr.segs order
%domain.showsegments(s, ones(size(s))); axis equal;

sq = domain(s([1 2 5 6 9 10 13 14]), 1); sq.deletecorner([1 3 5 7]);
si = [4 3 2 1 -1]; pm = [-1 -1 -1 -1 1];  % repeating indices for corner domains
ne = domain(s([4 3 2 1 15]), pm);
nw = domain(s(si+4), pm); sw = domain(s(si+8), pm); se = domain(s(si+12), pm);
e = domain([], [], s([16 12 8 4]), 1);      % exterior since starts w/ [] []
sq.setrefractiveindex(1.3);  % domain.showdomains([sq ne nw se sw e]);
s.setmatch('diel', 'tm');    % all segs get diel matching conds
pr = scattering([ne nw sw se e], [sq]);   % everything is 'air' except square

sq.addmfsbasis(@(t) R*exp(1i*t), 0.1, 90, []); % bit better than fb
for d = [ne nw sw se]          % attach nuFB to 4th corner each domain
  d.addnufbbasis(d.cloc(4), 2/3, d.cangoff(4), exp(-1i*pi/4)*d.cangoff(4), 15);
  d.addrpwbasis(20); %d.addregfbbasis(d.cloc(4), 15, []);
end
e.addmfsbasis(@(t) R*exp(1i*t), 0.2, 70, []); %pr.showbasesgeom; drawnow

pr.setoverallwavenumber(k);
pr.setincidentwave(pi/3);
tic; pr.solvecoeffs; fprintf('\tcoeffs done in %.2g sec\n', toc)
fprintf('\tL2 bdry error norm = %g, coeff norm = %g\n', ...
        pr.bcresidualnorm, norm(pr.co))
%pr.co = zeros(size(pr.co)); pr.co(200) = 1; % test individual basis funcs
figure; opts.dx = 0.05; opts.bb = [-2 2 -2 2];
tic; pr.showthreefields(opts); fprintf('\tgrid eval in %.2g sec\n', toc);
figure; imagesc(log10(abs(pr.A)));colormap(jet(256));caxis([-20 0]);colorbar;
title('A matrix for problem (log10 scale)');
