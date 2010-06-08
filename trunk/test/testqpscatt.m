% testing qpscatt class on various obstacles. Barnett 6/1/10
% Note: adapted from polesftlyp.m code (Dirichlet example)

clear all classes; v = 1;   % verbosity: 0 no pics, 1 final pic, 2 diagnostics
ob = 't';   % type of obstacle: 'd', 't', etc
N = 70; s = scale(segment.smoothstar(N, 0.3, 3), 0.25); % closed curve segment
om = 10; d = 1.0;              % incident wavenumber, problem x-periodicity
o.nei=0; o.buf=0; o.M = 150;       % M on FTy LPs (M>300 for Wood)
if ob=='d', de = domain([], [], s, -1);          % Dirichlet obstacle
  de.addlayerpot(s, 'd'); %[-1i*om 1]);          % adds CFIE to segment
  s.setbc(1, 'D', []);                           % homog Dirichlet BCs
  p = qpscatt(de, [], d, o);                     % create problem instance
elseif ob=='t', di = domain(s, 1); di.setrefractiveindex(1.5); % transmission
  de = domain([], [], s, -1);                    % exterior domain
  s.addinoutlayerpots('d');                      % new double-sided layerpots
  s.addinoutlayerpots('s');
  s.setmatch('diel', 'TM');                      % TM dielectric continuity
  p = qpscatt(de, di, d, o);  % airdoms must only be the single exterior domain
end
p.setoverallwavenumber(om);
p.setincidentwave(-pi/5); %-acos(1 - 2*pi/d.k)+1e-14; % 'single' Wood's anomaly
%for i=1:2, p.t.bas{i}.requadrature(p.t.bas{i}.N, struct('omega',om, 'nearsing', 3)); end  % how to enforce nearsing = 3; used to compare E blocks to std.mat
if v>1, p.showkyplane; end
p.fillquadwei;              % this is only for obstacle mismatch (blocks A,B)
% now maybe reset sqrtwei to 1 so vectors are plain values on bdry ?
p.fillrighthandside;
p.fillbcmatrix; fprintf('cond(E) = %.3g\n', cond(p.A))
if v>1, figure; imagesc(real(p.A)); set(gcf, 'name', 'testqpscatt: E matrix');
  colorbar; colormap(jet(256)); caxis(.01*[-1 1]);end

tic; p.solvecoeffs; toc;              % fills matrices, does least-squares soln
p.pointsolution(pointset(0.2+0.5i)) % ob='d': 0.090658571729 -0.014578950032i
p.bcresidualnorm
if v, p.showfullfield(struct('ymax', 1));title(sprintf('qpsc, obst=%s',ob));
  p.showbdry(struct('arrow',0,'normals',0,'label',0,'blobs',0)); end
if v>1, figure; p.showthreefields; h = findall(gcf,'Type', 'axes');
  set(gcf, 'CurrentAxes', min(h)); hold on;p.showbragg;end % adds Bragg to u_inc

