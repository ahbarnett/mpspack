% test the segment method invertZparam. Barnett 3/9/12

clear; N = 100; s = segment.smoothnonsym(N, 0.3, 0.3, 3);
d = domain(s, 1);  % for inside test

% test pts: grid of interior pts...
g = -1.3:0.02:1.3; [xx yy] = meshgrid(g); zg = xx(:)+1i*yy(:); % col, so z is
ii=d.inside(zg); z = zg(ii);

% line of pts fixed complex displacement from bdry...
%m = 1000; z = s.Z((1:m)/m + 0.4i/(2*pi));   % Im + is inside, - outside

% single point at known t value...
%to = (0 - 0.1i)/(2*pi); z = s.Z(to);

tic; t = s.invertZparam(z); toc  % do it (5 secs for grid dx=.05, n=100)

if size(t,2)==1, t*2*pi, end  % compare against number in to defn

figure; subplot(1,2,1); title('z plane (red dots should hit black ones)');
s.plot; hold on; plot(real(z),imag(z), 'ko');
plot(s.Z(t(find(~isnan(t)))), 'r.'); axis(1.5*[-1 1 -1 1]);
subplot(1,2,2); plot(2*pi*t, 'r.'); hold on; plot(2*pi*s.t, 0*s.t, '.-');
axis([0 2*pi -1 1]); %axis equal;
title('s = 2.pi.t plane');

if exist('xx') % do an interior contour plot of min imag part of s...
  t(find(imag(t)<0)) = nan+1i*nan; mis = nan*xx; % note nan for im crucial
  mis(ii) = min(imag(t), [], 1);
  figure; contour(g,g,mis); colorbar; hold on; s.plot;
end



% =======================
% independent test of domain.showimagparam:
clear; N = 100; s = segment.smoothnonsym(N, 0.3, 0.3, 3);
d = domain(s, 1);
o = []; o.levels = 0:0.02:0.12; [c h mis gx gy] = d.showimagparam(o);

% harder domain... with a bump, requires fine-tuning the o.to :
clear; N = 100;
a=.3;b=.3;w=3; a1=0.1; a2=0.1; t0=3*pi/4; % smoothnonsym w/ gaussian bump @ t0
s = segment.radialfunc(N, {@(q) 1 + a*cos(w*(q+b*cos(q))) - a1*exp(-0.5*((q-t0)/a2).^2), @(q) -a*sin(w*(q+b*cos(q))).*w.*(1-b*sin(q)) + (a1/a2^2)*(q-t0).*exp(-0.5*((q-t0)/a2).^2) });
d = domain(s, 1);
o = []; o.levels = 0:0.01:0.05; %o.dx = 0.1; % reduce grid density for speed
% following clusters start pts for invertZparam around gaussian bump...
n=10; no=10; o.to = [((1:n)-.5)/n, (t0+2*a2*(-no:2:no)/no)/(2*pi)] + 0.05i;
figure; [c h mis gx gy] = d.showimagparam(o);
figure; imagesc(gx,gy,mis); colorbar; set(gca,'ydir', 'normal');
hold on; s.plot;
