% testing qpscatt class on various obstacles, Bragg ampls. Barnett 6/1/10
% Note: adapted from polesftlyp.m code (Dirichlet example)

clear; v = 1;        % verbosity: 0 no pics, 1 final pic, 2 diagnostics
ob = 't'; d = 1.0;   % type of obstacle: 'd', 't', etc; problem x-periodicity 
N = 80; s = scale(segment.smoothstar(N, 0.3, 3), 0.35); % .25 closed curve
om = 10;             % incident wavenumber
o.nei = 2; o.buf = 1; o.M = 150;                 % M = # nodes per FTy LPs
if ob=='d', de = domain([], [], s, -1);          % Dirichlet obstacle
  de.addlayerpot(s, [-1i*om 1]);                 % adds CFIE to segment
  s.setbc(1, 'D', []);                           % homog Dirichlet BCs
  p = qpscatt(de, [], d, o);                     % create problem instance
elseif ob=='2', ss = [s translate(s.scale(0.6), .3-.6i)];
  de = domain([], [], num2cell(ss), {-1 -1});    % 2 obsts
  de.addlayerpot(ss, [-1i*om 1]);                 % adds CFIE to segment
  ss.setbc(1, 'D', []);                           % homog Dirichlet BCs
  p = qpscatt(de, [], d, o);                     % create problem instance
elseif ob=='t', di = domain(s, 1); di.setrefractiveindex(1.5); % transmission
  de = domain([], [], s, -1);                    % exterior domain
  s.addinoutlayerpots('d');                      % new double-sided layerpots
  s.addinoutlayerpots('s');
  s.setmatch('diel', 'TM');                      % TM dielectric continuity
  p = qpscatt(de, di, d, o);  % airdoms must only be the single exterior domain
end
%de.bas{1}.quad = 'a'; de.bas{2}.quad = 'a'; % tweak the quadrature
p.setoverallwavenumber(om);
p.setincidentwave(-pi/5);%-acos(1 - 2*pi/om));%-1e-14);%-pi/5;% 'single' Wood's anomaly
%for i=1:2, p.t.bas{i}.requadrature(o.M, struct('omega',om, 'nearsing', 3)); end  % how to enforce nearsing = 3; used to compare E blocks to std.mat
if v>1, p.showkyplane; end
p.fillquadwei;             % this is only for obstacle mismatch (blocks A,B)
p.sqrtwei = 1+0*p.sqrtwei; % row scaling: make vectors are plain values on bdry
p.fillrighthandside;
p.fillbcmatrix; fprintf('cond(E) = %.3g\n', cond(p.A))
if v>1, figure; imagesc(real(p.A)); set(gcf, 'name', 'testqpscatt: E matrix');
  colorbar; colormap(jet(256)); caxis(.01*[-1 1]);end

tic; p.co = p.linsolve; toc;    % least-squares soln
z = pointset(0.2+0.5i);
tic; p.pointsolution(z), toc % ob='d': 0.090658571729 -0.014578950028i
%A = p.evalbases(z); A*p.co % test if agrees with pointsolution
%[bo to bi ti] = p.braggampl(find(imag(p.kpy)==0));
%profile clear; profile on;
[u d n] = p.braggpowerfracs(struct('table',1)); %profile off; profile viewer

if v, tic; p.showfullfield(struct('ymax', 1)); toc
  title(sprintf('qpsc, obst=%s',ob));
  p.showbdry(struct('arrow',0,'normals',0,'label',0,'blobs',0)); end
if v>1, figure; p.showthreefields; h = findall(gcf,'Type', 'axes');
  set(gcf, 'CurrentAxes', min(h)); hold on;p.showbragg;end % adds Bragg to u_inc

%S = load('std'); figure; imagesc(abs(p.A - S.E)); % check E matrix correct
%figure; plot(eig(p.A), '+'); axis equal % show eigvals of E
