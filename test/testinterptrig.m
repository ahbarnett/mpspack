% test trig interpolation for upsampling by FFT. Barnett 3/9/12
% The key convention is the 0.5-gridpoint offset of both grids.

clear; n = 50; t = 2*pi*((1:n)-0.5)/n; % original grid
f = @(t) exp(2*sin(t)); % real analytic periodic
N = 3*n;
g = quadr.interptrig(f(t), N);
tn = 2*pi*((1:N)-0.5)/N;  % new grid
norm(g-f(tn))
figure; plot(t, f(t), '+'); hold on; plot(tn, g, '-');

% now test dense matrix for same thing:
A = quadr.interptrig(eye(n), N);
g = A * f(t).';
norm(g-f(tn).')

