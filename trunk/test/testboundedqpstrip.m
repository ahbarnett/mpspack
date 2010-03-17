% test bounded strip object, and qpstrip semi-infinite,
% ...and use them to qp a single source.
% barnett 3/6/10

clear all classes; v = 1;        % verbosity
thi = -pi/5; om = 10; %7.76644415490187; %(t.a=1) % inc ang, overall freq
%thi = -pi/5; om = 11.6496662323528; %(t.a=-1)
%thi = -pi/3; om = 4*pi; % 'double' Woods anomaly, with Bloch a=1
%om = 10; thi = -acos(1 - 2*pi/om);%'single' Wood's anomaly (generic Bloch a)
om = 2*pi/0.23; thi = 5*pi/8; % Yasumoto+Yoshitomi 1999 Table IV test case

%N = 35; M = 35; Mt = 15; L = 6; % Mt=top&bot # quadr, M=sides, L=Rayleigh order
N = 60; M = 46; Mt = 25; L = 8; % for Yasumoto99 test case
d = 1; yB = -1; yT = 1;       % strip & box geom
B = segment(Mt, 1i*yB+[-d/2,d/2], 'g'); T = segment(Mt, 1i*yT+[d/2,-d/2], 'g');
t = boundedqpstrip([B T], om, struct('M', M));
kvec = om*exp(1i*thi); 
a = exp(1i*real(conj(kvec) * d)); t.setbloch(a);
t.addregfbbasis(0.1+.1i, N, struct('rescale_rad', diam(t)));
t.setupbasisdofs; J = t.bas{1};  % J = Jexp basis
%figure; t.plot; t.showbasesgeom;
tt = qpstrip(d, om, struct('seg',T, 'pm',-1)); % create two semi-bounded strips
tb = qpstrip(d, om, struct('seg',B, 'pm',-1)); tt.setbloch(a); tb.setbloch(a);
tt.addqprayleighbasis(L); tb.addqprayleighbasis(L);
%figure;tt.plot(struct('gridinside',0.03));tb.plot(struct('gridinside',0.03));

tic; Q = t.evalbasesdiscrep; % fill matrix blocks and assemble system...
[Ac Anc] = tt.bas{1}.eval(T); [Ad And] = tb.bas{1}.eval(B);
W = [[Ac;Anc] zeros(2*Mt,2*L+1) ; zeros(2*Mt,2*L+1) [Ad;And]];
[AT ATn] = J.eval(T); [AB ABn] = J.eval(B); 
V = [AT;ATn;AB;ABn];
E = [-W V; zeros(2*M,2*(2*L+1)) Q];

nei = 2; % ......right-hand dide data: point source(s) nei=2 beats nei=1......
so = -.2+.3i; mfs = mfsbasis(pointset(so+(-nei:nei)'*d)); % sum of direct neigh
mfs.doms = t; % hook to bounded strip, now compute inhomog data and sum w/ phase
[AT ATn] = mfs.eval(T); [AB ABn] = mfs.eval(B); ph = a.^(-nei:nei).'; % phases
rhs = -[AT;ATn;AB;ABn;mfs.evalboundedstripdiscrep(t)]*ph;   % RHS
co = E \ rhs; fprintf('solution time %.3g s\n', toc);
fprintf('||E||_1 = %.3g, cond(E) = %.3g\n', norm(E,1), cond(E));
r = rhs - E*co; fprintf('norm resid for discrep = %.3g\n', norm(r));
coa = co(end-J.Nf+1:end); cocd = co(1:end-J.Nf); % pull out J, Rayleigh coeffs
fprintf('a coeff norm = %.3g, c&d coeff norm = %.3g\n', norm(coa), norm(cocd));
loc = so -0.2+0.0003i; % test Yasumoto1999 Table IV top block row (NB -x,y why?)
fprintf('u (Jexp+u_i) at (%g,%g) = %.16g\n',real(loc), imag(loc), ...
        t.evalbases(pointset(loc))*coa + mfs.eval(pointset(loc))*ph);
% cf Re[u]: y+0.03: 0.117120006144941, y+.003: 0.115891895634577
%           y+3e-4: 0.115881138140457. agree to 1e-15 error

if v, sc = 0.3; dx = 0.025; x = -1.5:dx:1.5; y = -2:dx:2;  % plot caxis, grid
[xx yy] = meshgrid(x,y); p = pointset(xx(:)+1i*yy(:));
uig = reshape(mfs.eval(p)*ph, size(xx));
figure; set(gcf,'position',get(gcf,'position') + [0 0 300 0]); subplot(1,3,1);
imagesc(x, y, real(uig)); hold on; t.plot; colormap(jet(256));
set(gca,'ydir','normal'); axis equal tight; caxis(sc*[-1 1]);
usg = reshape(t.evalbases(p)*co(end-J.Nf+1:end), size(xx));
title('u_{inc} to periodize'); subplot(1,3,2); imagesc(x, y, real(usg));
hold on; t.plot; set(gca,'ydir','normal'); axis equal tight; caxis(sc*[-1 1]);
title('QP field (Jexp)'); subplot(1,3,3); imagesc(x, y, real(uig+usg));
hold on; t.plot; set(gca,'ydir','normal'); title('total field');
axis equal tight; caxis(sc*[-1 1]);

figure; semilogy(abs(co), '+'); xlabel('j'); ylabel('|c_j|'); title('coeffs');
%figure; imagesc(real(W)); colorbar; caxis([-10 10])

end
