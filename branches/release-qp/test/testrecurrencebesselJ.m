% test accuracy of my fast Bessel code. Barnett Feb 2008.
% brought into mpspack, 7/23/09

clear all classes
if 1 % test small argument problems .............
x = [0 1e-18 1e-16 1.1e-16 1e-15 1e-10 1e-5 1 200]';   % note 100 makes nst>100
M = 2; %M = max(x);
Jex = besselJ(0:M, x);
J = utils.recurrencebesselJ(M, x);
J-Jex
end

% test speed ...............
Ms = [5 10 20 50 100 200];
for M = Ms;
disp(sprintf('M = %d .........................', M));
  x = M*sort(rand(1e4,1));    % as in app, use arguments out to largest order
  tic; Jex = besselJ(0:M, x); tex = toc;
  tic; J = utils.recurrencebesselJ(M, x); t = toc;
  nmev = (M+1)*numel(x)/1e6;  % how many time 1e6 bessel evals were done
  disp(sprintf('time (us per eval), matlab:%.3g, me:%.3g', tex/nmev, t/nmev));
  disp(sprintf('speedup factor %g\nmax error %g',tex/t,max(max(abs(J-Jex)))));
end
% I believe these errors are limited by matlab's besselj not my recurrence.

figure; imagesc(x,0:M,log10(abs(J-Jex)).');  xlabel('x'); ylabel('n');
caxis([-17 -14]); colormap(jet(256)); colorbar; set(gca,'ydir','normal');
title('abs err between matlab besselj and my fast one, on -log10 scale');


if 0    % .......... detailed timing, profiler
  M = 200; x = M*sort(rand(1e4,1));
  profile clear; profile on; J = utils.recurrencebesselJ(M, x); profile off;
  profile viewer;
end

if 0, figure;
  subplot(1,2,1); imagesc(x, n, Jex.'); xlabel('x'); ylabel('n');
  caxis([-.1 .1]); set(gca,'ydir','normal');
  subplot(1,2,2); imagesc(x, n, J.'); xlabel('x'); ylabel('n');
  caxis([-.1 .1]); set(gca,'ydir','normal');
end

