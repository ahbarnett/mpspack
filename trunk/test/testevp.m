% basic test routine for EVP class in MPSpack. Barnett 8/13/10
% Finds lowest 16 Dirichlet modes of smooth domain to 13 digits in 13 secs...

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

figure; imagesc(cumsum(s.w), 1:numel(p.kj), real(p.ndj)'); % plot bdry funcs
colormap(jet(256)); caxis([-1 1]*max(abs(caxis))); colorbar;
xlabel('s'); ylabel('j'); title('boundary functions \partial_n \phi_j (s)');
figure; plot(cumsum(s.w), real(p.ndj), '+-');
xlabel('s'); ylabel('\partial_n \phi_j (s)');
%evp.weylcountcheck(p.kwin(1), p.kj, d.perim, d.area);  % check missing?
