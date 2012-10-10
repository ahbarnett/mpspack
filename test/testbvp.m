% BVP class test routines
% barnett 7/18/08, included exterior domain and Neumann 7/19/08
% added MFS and transmission 7/20/08, added DLP 7/31/08, new interface 7/27/09

clear all classes
verb = 1;          % verbosity: 0 for no figures, 1 for figures
passdata = 0;      % 0 passes func handles, 1 passes data vectors
k = 10;            % overall wavenumber

for prob=1:9 % ===== loop over test problems ======
  
  switch prob      % choose segments, domains, basis sets, and exact solutions
   case {1,2,3,9} %    ...... interior polygon D & N (3: DLP, 9: test Greens th)
    dname = 'interior';
    if prob==3                            % choose a basis set for domain
      s = segment.smoothstar(70, 0.3, 3); pm = 1; d = domain(s, pm);
      d.addlayerpot([], 'd'); dname = [dname ' DLP'];
    elseif prob==9
      s = segment.smoothstar(70, 0.3, 3); pm = 1; d = domain(s, pm);
      d.addlayerpot([], 's'); d.addlayerpot([], 'd');
      dname = [dname ' Green Rep Thm'];
    else
      p = [1 1i exp(4i*pi/3)];            % triangle from Alex's inclusion paper
      M = 40;
      s(1) = segment(M, [p(1) p(2)]);
      s(2) = segment(M, [p(2) p(3)], 't');
      s(3) = segment(M, [p(1) p(3)]);       % note 3rd seg is backwards
      pm = [1 1 -1]; d = domain(s, pm);     % build interior polygonal domain
      opts.real = 1; d.addrpwbasis(30, opts);
    end
    xsing = 1 + 1i;                        % location of singularity
    u = @(x) bessely(0, k*abs(x - xsing)); % an exact Helmholtz soln in domain
    R = @(x) abs(x - xsing);
    du = @(x) -k*(x-xsing)./R(x).*bessely(1, k*R(x)); % grad u as C# (why? u=Re)
    ux = @(x) real(du(x)); uy = @(x) imag(du(x)); % u_x and u_y funcs
   
   case 4 % cases  ........... exterior polygon D & N Herglotz (ill-posed BVP)
    %segment.disconnect(s); s = s(end:-1:1);      % flip segment order for ext
    s(3) = segment(M, [p(1) p(2)]);
    s(2) = segment(M, [p(2) p(3)], 't');
    s(1) = segment(M, [p(1) p(3)]);
    pm = [1 -1 -1]; d = domain([], [], s, pm);   % new exterior polygonal domain
    dname = 'exterior Herglotz';
    N = 20; opts.real = 1; d.addrpwbasis(N, opts);
    u = @(x) besselj(0, k*abs(x));               % entire Helmholtz solution
    du = @(x) -k*x./abs(x).*besselj(1,k*abs(x)); % grad u as C#
    ux = @(x) real(du(x)); uy = @(x) imag(du(x)); % u_x and u_y funcs
   
   case {6,7} % ......................... ext analytic D & N radiative
    a = 0.3; w = 3;                             % use wiggly curve
    M = 80;
    R = @(t) 1 + a*cos(w*t); Rt = @(t) -w*a*sin(w*t);
    Z = @(s) exp(2i*pi*s).*R(2*pi*s);
    s = segment(M, {Z, @(s) 2*pi*(1i*Z(s) + exp(2i*pi*s).*Rt(2*pi*s))}, 'p');
    pm = -1; d = domain([], [], s, pm);
    dname = 'exterior radiative';
    N = 50; d.addmfsbasis(s, N, struct('tau',0.05,'real',1));
    xsing = -0.2 + 0.3i;                        % location of singularity
    u = @(x) bessely(0, k*abs(x - xsing)); % radiative Helmholtz soln
    R = @(x) abs(x - xsing);
    du = @(x) -k*(x-xsing)./R(x).*bessely(1, k*R(x)); % grad u as C# (why? u=Re)
    ux = @(x) real(du(x)); uy = @(x) imag(du(x)); % u_x and u_y funcs  
   
   case 8 % .............................. transmission, analytic curve
    a = 0.3; w = 3; M = 100;
    R = @(t) 1 + a*cos(w*t); Rt = @(t) -w*a*sin(w*t);
    Z = @(s) exp(2i*pi*s).*R(2*pi*s);
    s = segment(M, {Z, @(s) 2*pi*(1i*Z(s) + exp(2i*pi*s).*Rt(2*pi*s))}, 'p');
    d = domain(s, 1); d(2) = domain([], [], s, -1); % int, ext domains
    N = 50; d(2).addmfsbasis(s, N, struct('tau',0.05,'real',1)); opts.real=1;
    d(1).refr_ind=2; Ni = 40; d(1).addregfbbasis(0, Ni, opts);   % int basis set
    dname = 'transmission';
    t = 1; kvec = d(1).refr_ind*k*exp(1i*t);      % t is plane wave angle
    ui = @(x) exp(1i*real(conj(kvec) .* x)); % field inside (we take 0 outside)
    uix = @(x) 1i*real(kvec)*ui(x); uiy = @(x) 1i*imag(kvec)*ui(x); u = ui;
  end
  
  switch prob      % choose BCs
   case {1,3,4,6,9} % ................... Dirichlet
    bname = 'Dirichlet';
    for j=1:numel(s)                    % same func on all segments
      f = @(t) u(s(j).Z(t));            % compose u(Z(t))
      if passdata, f = f(s(j).t); end   % pass data vec via f instead of func?
      s(j).setbc(-pm(j), 'd', [], f);   % note -pm: which side domain BC is on
    end
   case {2,5,7} % ..................... Neumann
    bname = 'Neumann';
    for j=1:numel(s)                    % g will be normal deriv as func of t
      g = @(t) ux(s(j).Z(t)).*real(s(j).Zn(t)) + uy(s(j).Z(t)).*imag(s(j).Zn(t));
      if passdata, g = g(s(j).t); end   % pass data instead of func?
      s(j).setbc(-pm(j), 'n', [], g);
    end
   case 8 % ............................ Matching condition
    bname = 'Matching';
    f = @(t) -ui(s.Z(t));             % jump in value is 0-ui (outside-inside)
    g = @(t) -uix(s.Z(t)).*real(s.Zn(t)) - uiy(s.Z(t)).*imag(s.Zn(t)); % jump un
    if passdata, f = f(s.t); g = g(s.t); end
    s.setmatch([1 -1], [1 -1], f, g);    % [u]=[un]=0 matching
  end

  pr = bvp(d); % ......... set up then solve BVP, plot soln and err
  pr.setoverallwavenumber(k);
  tic; if prob==9 % GRF test, no lin solve, just set up RHS by hand
      un = real(conj(du(s.x)).*s.nx); pr.co = [un; -u(s.x)];
      pr.fillquadwei; pr.fillbcmatrix; pr.rhs = u(s.x).*pr.sqrtwei.';
     else pr.solvecoeffs; end
  name = sprintf('BVP #%d, %s %s', prob, dname, bname);
  fprintf('%s: done in %.2g sec; cond(A)=%2.g\n', name, toc, cond(pr.A))
  fprintf('\tL2 bdry error norm = %g, coeff norm = %g\n', pr.bcresidualnorm, ...
          norm(pr.co))
  if verb, figure('name', name); subplot(1,2,1);
  opts.dx = 0.03; tic; [uN gx gy di] = pr.gridsolution(opts); t = toc;
  imagesc(gx, gy, real(uN)); title('bdry and soln uN');
  c=caxis; caxis(max(c)*[-1 1]);
  axis equal tight; colorbar; set(gca,'ydir','normal'); hold on; pr.showbdry;
  [xx yy] = meshgrid(gx, gy); zz = xx + 1i*yy;
  ud = u(zz); ii = find(~isnan(di)); % ii = indices which are inside any domain
  if prob==8, ud(find(di==2)) = 0; end % transm: kill exact value outside
  subplot(1,2,2); imagesc(gx, gy, real(uN - ud)); c=caxis; caxis(max(c)*[-1 1]);
  axis equal tight; colorbar;
  set(gca,'ydir','normal'); title('soln err: uN - u');
  r = opts.dx * norm(uN(ii) - ud(ii));    % L2 domains error norm (est on grid)
  fprintf('\tgrid eval in %.2g sec, L2 domain error norm = %g\n\n', t, r);
  end
end % ====== end prob loop ========
disp('All bdry err norms should be small, and all domain norms small except #3 and #9. These errs are large due to the fact that no special scheme was used for close evaluation of layer-pots.');
