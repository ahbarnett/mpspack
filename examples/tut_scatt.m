% Example codes from MPSpack tutorial, also generates EPS figures for this doc
% SECTION 6: SCATTERING

clear all classes; verb = 0;          % if verb>0, generates EPS figures

tref = segment.radialfunc(250, {@(q) 1 + 0.3*cos(3*q), @(q) -0.9*sin(3*q)});
d = domain([], [], tref, -1);                    % overwrites previous d
tref.setbc(1, 'D', []);       % homogeneous Dirichlet BCs: sound-soft
opts.tau = 0.05; d.addmfsbasis(tref, 200, opts);  % basis set for ext domain
p = scattering(d, []);
k=30; p.setoverallwavenumber(k);
p.setincidentwave(pi/6);              % incident plane wave with given angle
tic; p.solvecoeffs; toc;              % fills matrices, does least-squares soln
p.bcresidualnorm
p.showthreefields;

if verb   %   generate f:soft
  figure; o.dx=0.01; o.bb = 2*[-1 1 -1 1];
  p.showthreefields(o); %set(gca,'fontsize',20);
  print -depsc2 ../doc/soft.eps
end

tref.setbc(1, 'N', []);       % homogeneous Neumann: sound-hard
tref.setbc(1, 1i*k, 1);       % homogeneous Robin: impedance

% transmission problem
di = domain(tref, 1); di.refr_ind = 1.5;
opts.tau = -0.03; di.addmfsbasis(tref, 220, opts);
tref.setmatch('diel', 'TM');
p = scattering(d, di); p.setoverallwavenumber(k); p.setincidentwave(pi/6);
%figure; p.plot;
p.solvecoeffs; p.bcresidualnorm, p.showthreefields;

if verb   %   generate f:diel
  figure; o.dx=0.01; o.bb = 2*[-1 1 -1 1];
  p.showthreefields(o); %set(gca,'fontsize',20);
  print -depsc2 ../doc/diel.eps
end
