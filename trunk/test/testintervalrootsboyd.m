% test Boyd spectral rootfinding
% Barnett 7/28/2010, 8/13/10

f = @(x) sin(pi*x);

disp('high accuracy request...')
[x e y u ier] = utils.intervalrootsboyd(f, [0 10]); x, max(e)
figure; plot(y, u, '-'); hold on; plot(x, 0*x, 'rx'); title('intervalrootsboyd')

disp('low accuracy request...')
[x e y u ier] = utils.intervalrootsboyd(f, [0 10], struct('Ftol',1e-6)); x, max(e)
size(u)

disp('taxing one...')
[x e y u ier] = utils.intervalrootsboyd(f, [0 100]); x, max(e)
% note that errors e are too pessimistic near endpoints of interval - why?

disp('should fail gracefully...')
[x e y u ier] = utils.intervalrootsboyd(f, [0 1e3]); ier
