% test insertion of basis sets within domains, and plotting inside domains
% barnett 7/8/08

clear all classes
M = 20; s = segment.polyseglist(M, [1 1i exp(4i*pi/3)]);  % CCW tri from INCL
d = domain(s, 1);
k = 10; opts.real = 1; opts.fast = 0;
d = d.addregfbbasis(0, 10, k, opts); % must return d since domain is value class
js = 1:d.Nf;                       % which basis func indices to plot

[zz ii gx gy] = d.grid(0.01);      % set up grid then evaluate
A = d.evalbases(pointset(zz));

n = numel(js);                     % ... now show the basis sets 
nh = floor(sqrt(1.8*n)); nv = ceil(n/nh); % how many across and down, subplot
figure('name', 'reg FB basis: Re[u_j(x)]'); c = .5;    % caxis range
for i=1:numel(js)
  subplot(nv, nh, i); u = NaN*ones(numel(gy), numel(gx)); % prepare plot grid
  u(ii) = A(:,js(i)); imagesc(gx, gy, real(u));
  caxis(c*[-1 1]); set(gca, 'ydir', 'normal'); axis off equal tight;
end
subplotspace('vertical', -30); subplotspace('horizontal', -30);
d.plot;
