% low eigenvalues of non-simply connected domain via Fredholm determinant method
% Barnett 2011, w/ Braxton Osting.

clear % one analytic segment, periodic radial function
a = 0.3; w = 3;
M = 240;
R = @(t) 1 + a*cos(w*t); Rt = @(t) -w*a*sin(w*t); Rtt = @(t) -w^2*a*cos(w*t);
s = [];
s{1} = segment.radialfunc(M, {R, Rt, Rtt});  % easier than you
s{2} = segment(30, [.5+0.0i, .3, 0, 2*pi], 'p');  % note I added to s not si
d = domain(s{1},1, s{2}, -1);
%figure (1); d.plot; drawnow;

% note +-1 is wrt segment's natural sense, not their sense with domain list:
s{1}.setbc(-1, 'D'); s{2}.setbc(1, 'D');

% note if you don't supply Rtt above, Alpert quad has to be used not Kress:
for ii=1:numel(s); d.addlayerpot(s{ii}, 'd'); d.bas{ii}.quad = 'm'; end

p = evp(d);
opts = struct('modes',1,'verb',1,'tol',1e-12,'iter',1);
tic; p.solvespectrum([2 12], 'fd', opts); toc
figure; p.weylcountcheck(p.kj(1),p.kj,d.perim,d.area);
figure; p.showmodes
% problem is degeneracies only one mode per eigenspace computed.
