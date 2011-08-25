% Test and demo routine for EVP solvespectrum 'ms' in MPSpack. Barnett 8/23/11

clear all classes;
s = segment.smoothnonsym(130, 0.3, 0.2, 3); % create a closed segment
d = domain(s, 1);                           % create an interior domain
s.setbc(-1, 'D');                           % put Dirichlet BCs on inside
p = evp(d);                                 % create eigenvalue problem
d.addlayerpot(s, 'd');  % needed for 'fd' and 'ms' but not 'ntd'
kint = [2.5 9]; p.solvespectrum(kint, 'fd'); kjgood = p.kj;  % reference method

if 0, ks = 100:0.005:100.1; % check slopes of sing vals of (1/2-D): max is 1.5
  for i=1:numel(ks), ss(:,i)=svd(p.fillfredholmop(ks(i))); ss(end,i), end
  figure; plot(ks, ss, '+-'); end

o.maxslope = @(k) 1.5; o.verb = 0;
tic; p.solvespectrum(kint, 'ms', o); toc, fe = p.err.mininfo.fevals;
fprintf('fevals = %d  (%g per min found)\n', fe, fe/numel(p.kj)) 
p.kj - kjgood
p.err.ej


% to use...

%tic; p.solvespectrum([2.5 9], [], struct('modes',1)); toc  % solve everything

figure; imagesc(cumsum(s.w), 1:numel(p.kj), real(p.ndj)'); % image bdry funcs
colormap(jet(256)); caxis([-1 1]*max(abs(caxis))); colorbar;
xlabel('s'); ylabel('j'); title('boundary functions \partial_n \phi_j (s)');
figure; plot(cumsum(s.w), real(p.ndj), '+-');  % plot them as overlayed graphs
xlabel('s'); ylabel('\partial_n \phi_j (s)');

if 0, evp.weylcountcheck(p.kwin(1), p.kj, d.perim, d.area); end % check missing?

tic; p.showmodes; toc                 % compute and plot all modes

if 0, p.showmodes(struct('inds',[1 16 2])); % test choosing modes, by index...
  p.showmodes(struct('kwin',[7 8]));  % ...and by wavenumber window
end

if 0, [uj gx gy di] = p.showmodes;  % check output and normalization (GRF eval)
  for j=1:numel(p.kj); u=uj(:,:,j); sum(u(find(di==1)).^2)*(gx(2)-gx(1))^2, end
end  % squared L2-norms shown should be within O(grid spacing) of unity

