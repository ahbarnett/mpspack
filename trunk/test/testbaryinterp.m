% test for @utils/bary*.m Lagrange barycentric intepolation
% barnett 1/27/11

sc = 1.0; %1e10;  % scale (causes underflow, etc...?)
x = (1:10)/sc;   % interp nodes
gx = (2:0.1:9)/sc;  % evaluation checking nodes
f = @(x) sin(x*sc);
w = utils.baryweights(x);
y = f(x); gy = f(gx);
%y = 0*y; y(5) = 1;  % extract 5th col of L
u = utils.baryeval(x, w, y, gx);
figure; plot(x, y, '+'); hold on; plot(gx, [gy; u], '-');
figure; plot(gx, gy - u, '-'); title('difference');

