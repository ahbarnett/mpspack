% test Boyd spectral rootfinding
% Barnett 7/28/2010

f = @(x) sin(pi*x);
[x e y u] = utils.intervalrootsboyd(f, [0 10]); x, max(e)
figure; plot(y, u, '-'); hold on; plot(x, 0*x, 'rx'); title('intervalrootsboyd')

[x e y u] = utils.intervalrootsboyd(f, [0 10], struct('Ftol',1e-6)); x, max(e)
size(u)
