% DIELHOLE: scattering from dielectric with 'metallic' and air inclusions.
% Barnett 7/27/09, taken from 7/27/08 test/testscattering.m case 7.

clear all classes;
s = segment.smoothstar(200, 0.2, 3);            % weak trefoil shape (radius~1)
c = segment.smoothstar(70, 0.1, 2);             % squashed circle (radius~1)

sd = s.scale(2);                                % outer boundary of dielectric
sm = c.scale(0.5);
sa = translate(rotate(sm,.3), .8);              % air inclusion boundary
sm.rotate(pi/5); sm.translate(-.8-.6i);         % 'metal' boundary

de = domain([], [], sd, -1);                    % exterior air domain
da = domain(sa, 1);                             % air inclusion domain
d = domain(sd, 1, {sm sa}, {-1 -1});            % dielectric w/ 2 holes
d.setrefractiveindex(1.5);                      % choose dielectric index

setmatch([sd sa], 'diel', 'TE');                % impose TE matching diel-air
sm.setbc(1, 'D', []);                           % Dirichlet (PMC) BC on 'metal'

pr = scattering([de da], d);          % includes air pocket in nonzero u_inc
pr.setoverallwavenumber(8);
pr.setincidentwave(-pi/3);    % if just angle given, assumed plane wave
 
de.addmfsbasis(sd, 120, struct('tau',0.05)); % BASIS SETS: exterior domain
da.addmfsbasis(sa, 60, struct('tau',-0.1)); % inside air pocket
d.addmfsbasis(sd, 150, struct('tau',-0.05)); % 1st of 3 curves for diel...
d.addmfsbasis(sm, 60, struct('tau',0.1));   % ...2nd curve in metal pocket
d.addmfsbasis(sa, 60, struct('tau',0.1));   % ...3rd curve in air pocket

tic; pr.solvecoeffs; fprintf('\tcoeffs done in %.2g sec\n', toc)
pr.bcresidualnorm

pr.pointsolution(pointset(0))               % solution u value at origin

if 0          % output the figure used in the tutorial PDF file
  opts.dx = 0.02; opts.bb = [-4 4 -4 4]; pr.showthreefields(opts);
  print -depsc2 ../doc/figs/dielhole.eps
  figure; pr.plot; pr.showbasesgeom; axis tight;
  print -depsc2 ../doc/figs/dielholegeom.eps
end

% SET UP CONVERGENCE STUDY...

sd.requadrature(400); sm.requadrature(140); sa.requadrature(140);
d.clearbases; de.clearbases; da.clearbases;
de.addmfsbasis(sd, [], struct('tau',0.05,'nmultiplier', 120/450));
da.addmfsbasis(sa, [], struct('tau',-0.1,'nmultiplier', 60/450));
d.addmfsbasis(sd, [], struct('tau',-0.05, 'nmultiplier', 150/450));
d.addmfsbasis(sm, [], struct('tau',0.1, 'nmultiplier', 60/450));
d.addmfsbasis(sa, [], struct('tau',0.1, 'nmultiplier', 60/450));
pr.updateN(450);
tic; pr.solvecoeffs; fprintf('\tcoeffs done in %.2g sec\n', toc)
pr.bcresidualnorm

Ns=300:50:600; for i=1:numel(Ns)
  pr.updateN(Ns(i)); pr.solvecoeffs; u0 = pr.pointsolution(pointset(0));
  fprintf('N=%d, co-norm=%g, bc-norm=%g, |u(0)|=%.16g\n', Ns(i), ...
          norm(pr.co), pr.bcresidualnorm, abs(u0));
end
