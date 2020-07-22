% Helmholtz scattering through many lenses. For ommatidium eye project.
% (Ali Simons, Leslie Greengard, Charlie Epstein.)
% Alex Barnett 7/22/20
clear; clear all classes; addpath ..
verb = 1;
shape = 'connectdots';

switch shape
 case 'ellipse'         %    analytic ellipses, semiaxes a,b
  a = 0.3; b = 1;
  pp = 2*pi;
  Z = @(t) a*cos(pp*t)+1i*b*sin(pp*t);
  Zp = @(t) pp*(-a*sin(pp*t)+1i*b*cos(pp*t));
  Zpp = @(t) pp^2*(-a*cos(pp*t)-1i*b*sin(pp*t));
  n = 160;     % discr nodes for 10 digits (check conv of soln at test pt!)
  s0 = segment(n,{Z,Zp,Zpp},'p');
 case 'connectdots'     % arbitrary, connect the dots
  np = 40;    % # dots   (shoudl be multiple of 4)
  z = exp(2i*pi*(-np/4:np/4-1)/(4*np)) - cos(pi/8);     % a pi/8 angle lens
  z = [z -z];  % some points, replace w/ yours
  z = 3*z;
  n = 160;     % discr nodes for 5 digits.
  s0 = segment.smoothfourierz(n,z,n/2,1e-4);   % see help for this command
 case 'twoarcs'
  % *** this would have corners. I can easily do this if you need.
end
  
% make 3 lenses... (cell array since separate objects)
s{1} = s0;
s{2} = translate(s0,1);
s{3} = translate(s0.scale(0.5),2);
% interior domains, refr indices...
di(1) = domain(s{1}, 1); di(1).setrefractiveindex(1.3);
di(2) = domain(s{2}, 1); di(2).setrefractiveindex(1.5);
di(3) = domain(s{3}, 1); di(3).setrefractiveindex(2.0);
de = domain([], [], s, {-1 -1 -1});  % exterior (cell arr -> discon)
if verb, figure; di.plot; hold on; de.plot; axis equal; drawnow; title('geom');
end

% rep...
o.quad = 'm';                                     % Kress spectral quadr
for k=1:numel(s)           % add Kress-Roach rep...
  s{k}.addinoutlayerpots('d', o);      % rep affects both sides.
  s{k}.addinoutlayerpots('s', o);      % "
  setmatch(s{k}, 'diel', 'TM');        % matching conds: [u]=0, [u_n]=0
end
pr = scattering(de, di);                          % set up scattering BVP

k=20; pr.setoverallwavenumber(k);
pr.setincidentwave(-0.05*pi);     % if just angle given, it's a plane wave
% compute...
pr.fillquadwei; pr.setupbasisdofs;
pr.fillrighthandside;
pr.fillbcmatrix;
pr.linsolve;

ptest = pointset(2.5+0.1i);
pr.pointsolution(ptest)        % check u_scatt at one exterior pt

if verb, figure; opts.dx = 0.02; opts.bb = [-1 3 -2 2];
tic; pr.showthreefields(opts); fprintf('\tgrid eval in %.2g sec\n', toc);
hold on; plot(ptest.x,'.','markersize',20);    % test pt
end
