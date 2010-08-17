% Test and demo routine for EVP class in MPSpack. Barnett 8/17/10
% Finds lowest 16 Dirichlet modes of smooth domain to 13 digits in 13 secs...
% then plots them on a grid. See end for demo of finding more modes.

clear all classes; a = 0.3; b = 0.2; w = 3;  % shape params, smooth closed curve
s = segment.radialfunc(160, {@(q) 1 + a*cos(w*(q+b*cos(q))), ...
                    @(q) -a*sin(w*(q+b*cos(q))).*w.*(1-b*sin(q)), ...
                    @(q) -a*cos(w*(q+b*cos(q))).*w^2.*(1-b*sin(q)).^2 + ...
                    a*sin(w*(q+b*cos(q))).*w.*b.*cos(q)});  % includes curvature
d = domain(s,1);
s.setbc(-1, 'D');             % homog dirichlet BCs (on inside: note -1)
d.addlayerpot(s, 'd');
p = evp(d);                   % sets up problem object
tic; p.solvespectrum([2.5 9], [], struct('modes',1)); toc  % solve everything

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


% More impressive demo: same shape, modes 1-93 found and plotted in 100 seconds
tic; p.solvespectrum([2.5 20], [], struct('modes',1)); toc      
tic; p.showmodes; toc
max(p.err.minsigj) % since N=160 still, see some deterioration approaching k=20
