% test routine for layer potential basis objects
% Also see testlpquad for more low-level tests of layerpot.S etc.
% barnett 7/30/08

clear all classes
sh = 't';                    % 't' or 's' for triangle or smooth curve

if sh=='s', s = segment.smoothstar([], 0.3, 3);  % periodic analytic
elseif sh=='t', s = segment.polyseglist([], [1 1i exp(4i*pi/3)]);  % triangle
end
k = 0;                     % allows tau=1 test to work
d = domain(s, 1);                    % interior domain
[g.x g.ii g.gx g.gy] = d.grid(0.05);                     % some interior points
Ms = 100:50:100
for m = 1:numel(Ms)
  M = Ms(m); s.requadrature(ceil(M/numel(s))); % can also try 'g', better
  numel(vertcat(s.w))          % show # quad pts total
  d.clearbases; d.addlayerpotbasis([], 'd', k); D = d.evalbases(g);
  tau = ones(size(D,2),1);      % unity col vec
  ug = D * tau;                 % field on interior grid pts
  figure; d.plot; axis equal; u = NaN*zeros(size(g.ii)); u(g.ii)=ug;
  imagesc(g.gx, g.gy, log10(abs(u+1))); colorbar; caxis([-16 0]);
end

% interior Helmholtz GRF via SLP+DLP............................................
k = 0;
d.addlayerpotbasis([], 'd', k);
o.layerpotside = -1;                 % the inside
D = [];
for j=1:numel(d.bas), D = [D d.bas{j}.eval(g, o)]; end
% check interior pts have value -1
tau = ones(size(D,2),1);      % unity col vec
ug = D * tau;                 % field on interior grid pts
% ... unfinished...

