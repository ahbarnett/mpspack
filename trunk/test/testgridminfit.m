% test the gridded minimum fit routine from @evp
% barnett 8/23/11. Use o.verb=3 for figure-plotting each recursion step.

clear; close all; o.verb = 1; o.xtol = 1e-12;
f = @(x) abs(sin(x));  % number of vector entries m=1 
g = 0:1:100;  % grid
[x y i] = evp.gridminfit(f, g, o);
x/pi, y    % should be integers and 0 for values
fprintf('fevals = %d  (%g per min found)\n', i.fevals, i.fevals/numel(x))

f = @(x) abs([sin(x); sin(x+pi/5)]);  % number of vector entries m=2
o.verb = 3; g = 0:1:20;  % grid
[x y i] = evp.gridminfit(f, g, o);
x/pi, y
fprintf('fevals = %d  (%g per min found)\n', i.fevals, i.fevals/numel(x))
hold on; xx=0:0.03:20; plot(xx,f(xx),'r-'); plot(x,y(1,:),'go');

f = @(x) abs([sin(x-pi/10); sin(x+pi/5)]); % now try recursive, double mins
o.verb = 1; o.maxslope = @(x) 1; g = 0:1:20; %0:1:20;  % grid
[x y i] = evp.gridminfit(f, g, o);
x/pi, y
fprintf('fevals = %d  (%g per min found)\n', i.fevals, i.fevals/numel(x))
figure(1); hold on; xx=min(g):0.03:max(g); plot(xx,f(xx),'r-');
plot(x,y(1,:),'go');

f = @(x) abs([sin(x); sin(x+pi/100)]); % closer double mins: works harder
o.verb = 0; o.maxslope = @(x) 1; g = 0:.5:20; %0:1:20;  % grid
[x y i] = evp.gridminfit(f, g, o);
x/pi, y
fprintf('fevals = %d  (%g per min found)\n', i.fevals, i.fevals/numel(x))
figure(1); if o.verb>=3, hold on; end
xx=min(g):0.03:max(g); plot(xx,f(xx),'r-'); plot(x,y(1,:),'go');
% 27 fevals per min

% best model for min-svd ..................................................
N=100; b=sort(rand(1,N)); f = @(x) abs(exp(x-b)-1); % curved V-shapes in [0,1]
o.verb = 0; o.maxslope = @(x) 1.0; g = 0:0.2/N:1; % init grid (5 per typ gap)
[x y i] = evp.gridminfit(f, g, o);
if numel(x)==numel(b), max(abs(x-b)), else x, b, end  % g=(1:4);b(find(b<g(4)))
fprintf('fevals = %d  (%g per min found)\n', i.fevals, i.fevals/numel(x))
figure(1); if o.verb>=3, hold on; end
xx=0:1e-3:1; yy=nan(N,numel(xx)); for i=1:numel(xx), yy(:,i)=f(xx(i)); end
plot(xx,yy,'r-'); hold on; plot(x,y(1,:),'go'); axis([0 1 0 0.2]);
% with maxslope correct, about 17 evals per min
% maxslope 1/2 it's true size: misses around 10% of mins
% maxslope 2x correct, 21 evals per min (init grid 5x mean spacing)
