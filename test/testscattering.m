% scattering class test routine, including dielectrics, metallic bdry, etc
% barnett 7/27/08
% NEEDS TO BE CHANGED TO NEW k-FREE BASIS INTERFACE, MFS.

clear classes
opts = []; verb = 1;  % verbosity: 0 for no figures, 1 for figures (slower)
k = 8;                % note M,N do not yet adapt to different k
M = 100; s = segment.smoothstar(M, 0.2, 3);  % a recurring closed segment

for prob=6  % ======= scenarios
  s.disconnect;              % lets us reuse segment s afresh
  opts.testtransparent = 0;
  switch prob
    
   case {1,2},                     name = 'single metallic'; % ..............
    de = domain([], [], s, -1);
    if prob==1, bc='D'; else bc='N'; end; name = [name ' ' bc];
    s.setbc(1, bc, []);
    pr = scattering(de, []);  % note basis sets can be added after problem set
    N = 70; bopts.real = 0;   % complex (radiative) MFS basis
    de.addmfsbasis([], 0.4, N, [], bopts);  % no k needed
   
   case {3,4,5},                   name = 'single dielectric'; % ...........
    de = domain([], [], s, -1);
    transp = prob>4;    % if true, show u_err. must be false for non-diel case
    opts.testtransparent = transp;
    if transp, n = 1.0; name = [name ', transparent']; else, n = 1.2; end
    d = domain(s, 1); d.setrefractiveindex(n);
    if prob==3, pol='TM'; else pol='TE'; end; name = [name ' ' pol];
    s.setmatch('diel', pol); pr = scattering(de, d);
    N = 60; de.addmfsbasis([], 0.5, N, []);
    if transp, d.addregfbbasis(0, 40, []);
    else                       % for diel use real MFS ins
    bopts.real = 0; d.addmfsbasis([], -0.4, 80, [], bopts); end
    
   case {6,7},                   name = 'diel w/ air & metal pockets'; % .....
    s = segment.smoothstar(200, 0.2, 3); % original w/ more M
    c = segment.smoothstar(70, 0.1, 2); % squashed circle
    sm = c.scale(0.5); sd = s.scale(2); sa = translate(rotate(sm,.3), .8);
    sm.rotate(pi/5); sm.translate(-.8-.6i);
    de = domain([], [], sd, -1); da = domain(sa, 1); % exterior, air pocket
    d = domain(sd, 1, {sm sa}, {-1 -1}); d.setrefractiveindex(1.2);
    if prob==6, pol='TM'; else pol='TE'; end; name = [name ' ' pol];
    setmatch([sd sa], 'diel', pol); sm.setbc(1, 'D', []);
    pr = scattering([de da], d);          % incl the middle bit (air) in u_inc
    de.addmfsbasis([], 0.2, 110, []);
    da.addmfsbasis([], -0.6, 40, []); %da.addregfbbasis(.8, 30); %for air pocket
    bopts.segnum = 1; d.addmfsbasis([], -0.3, 150, [], bopts); % ext + 2 holes
    bopts.segnum = 2; d.addmfsbasis([], 0.7, 50, [], bopts);
    bopts.segnum = 3; d.addmfsbasis([], 0.7, 50, [], bopts);
 end
  pr.setoverallwavenumber(k); % note may be done before or after bases chosen
  pr.setincidentwave(pi/2 - pi/20);  % if just angle given, it's a plane wave
  fprintf('Prob #%d: %s\n', prob, name)
  tic; pr.solvecoeffs; fprintf('\tcoeffs done in %.2g sec\n', toc)
  fprintf('\tL2 bdry error norm = %g, coeff norm = %g\n', ...
          pr.bcresidualnorm, norm(pr.co))
  opts.dx = 0.05; opts.bb = [-3 3 -3 3];
  if verb, tic;
      figure('name', sprintf('prob #%d: %s', prob, name), 'position', [400 400 900 300]);
      pr.showthreefields(opts); fprintf('\tgrid eval in %.2g sec\n', toc);
  end
end % =========

figure; pr.showbdry; pr.showbasesgeom; axis equal; title('bdry & MFS pts');

%if verb>1, pr.showbdry; domain.showdomains(d); axis equal; end
