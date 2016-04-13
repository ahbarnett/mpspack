% accessibility modes in Fig.1 w/ Tacy, Hassell. Barnett 10/28/15
% Uses Fred det method. Takes 1 min to run.
clear

if 1 % ------ disc
  M = 400;
  s = segment(M,[0 1.0 0 2*pi], 'p');
  d = domain(s, 1);        % create an interior domain
  p = evp(d);              % sets up eigenvalue problem object
  s.setbc(-1, 'N');
  d.addlayerpot(s, 's');          % SLP set appropriate for Neu BC
  n = 30;
  kn = fzero(@(x) besselj(n-1,x) - n*besselj(n,x)./x, n)   % root of Jn'
  %k = 32.5342235567901
  o.tol = 1e-7; % degenerate, so the tolerance needs to be loosened... (sadly)
  o.modes = 1; p.solvespectrum(kn+0.05*[-1 1], 'fd', o);
  o.inds = 1; o.dx = 0.003; p.showmodes(o);
  text(-1,1,'(a)'); set (gcf,'paperposition',[0 0 2.5 2.5])
  print -depsc2 disc_tataru.eps
end

if 1   %----- generic smooth
  M = 450; % 1000 to check converged
  a = 0.3; b = 0.2; w = 3; % shape params (default for BNDS: a=.3 b=.2 w=3)
  s = segment.radialfunc(M, {@(q) 1 + a*cos(w*(q+b*cos(q))), ...
                      @(q) -a*sin(w*(q+b*cos(q))).*w.*(1-b*sin(q)), ...
                      @(q) -a*cos(w*(q+b*cos(q))).*w^2.*(1-b*sin(q)).^2 + a*sin(w*(q+b*cos(q))).*w.*b.*cos(q)});  % includes curvature. 200
  d = domain(s, 1);        % create an interior domain
  p = evp(d);              % sets up eigenvalue problem object
  s.setbc(-1, 'N');
  d.addlayerpot(s, 's');          % SLP set appropriate for Neu BC
  k = 40;
  o.modes = 1; p.solvespectrum([40.5 40.52], 'fd', o);
  o = []; o.inds = 1; o.dx = 0.005; p.showmodes(o);
  %k=40.51282199500848   (M=450)
  %  40.51282199500847   (M=1e3)
  text(-1,1.2,'(b)'); set (gcf,'paperposition',[0 0 3 3])
  print -depsc2 neu_mode.eps
end

