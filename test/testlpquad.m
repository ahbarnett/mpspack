% test Kapur-Rokhlin, Kress, quadrature on single & double layer potentials
% barnett 2/5/08 adapted to MPSpack 7/31/08
clear classes
verb = 0;            % verbosity

% ================== test SLP/DLP vals, and in Laplace case const DLP Greens
k = 0;                      % wavenumber (>=0 for 's' or 'd', =0 for 'c')
ords = [2 6 10 Inf];         % Kapur orders (Inf gives Kress)
co = 'bgrk';                 % plot colors
test = 'c';                  % 'd'/'s' for DLP/SLP vals, 'c' for const DLP k=0

if test=='c', figure; end
for o = 1:numel(ords) % -------------------------- loop over orders
  ord = ords(o)
  opts.quad = 'k'; opts.ord = ord; if ord==Inf, opts.quad = 'm'; end
  
  Ns=50:50:450; err = zeros(size(Ns));        % quadrature convergence
  for i=1:numel(Ns);
    N = Ns(i);
    s = scale(segment.smoothstar(N, 0.3, 5), 0.5); % same pentafoil as in pc2d
    if test=='s', S = layerpot.S(k, s, [], opts);
    else, D = layerpot.D(k, s, [], opts); end
    tau = ones(N,1);              % const DLP density should give -1 inside
    if test=='c'
      err(i) = norm((D - eye(N)/2)*tau + 1)*sqrt(2*pi/N); %  rough bdry L2 norm
    elseif test=='d'
      u = (D - eye(N)/2)*tau;        % val, includes jump relation for inside
      disp(sprintf('k=%g N=%d\tDLP val at y_N = %.16g',k, N, u(end)));
    else
      u = S*tau;                     % val, no jump relation for val
      disp(sprintf('k=%g N=%d\tSLP val at y_N = %.16g',k, N, u(end)));
    end
  end
  if test=='c'
    loglog(Ns, err, ['+-' co(o)]);
    hold on; plot(Ns, 1e-2*(Ns/1e2).^-ord, ['--' co(o)]); drawnow % slope
  end
end % -----------------------------------------
if test=='c', axis tight; legend(num2cellstr(kron(ords, [1 1])));
  xlabel('N'); ylabel('L_2 err');
  title(['DLP k=0 tau\equiv1 convergence at orders: ' sprintf('%d ', ords)]);
  if verb>1, print -depsc2 test_dlpquad_laplace.eps
  end
end
if verb, figure; imagesc(real(D)); colorbar; figure; show_bdry(bd); end


% ==================== Test Greens Thm for combined SLP/DLP, k>=0
k = 10;                      % wavenumber (may be 0 or >0)
ords = [2 6 10 Inf];         % Kapur orders (Inf gives Kress)
co = 'bgrk';

figure;
for o = 1:numel(ords) % -------------------------- loop over orders
  ord = ords(o)
  opts.quad = 'k'; opts.ord = ord; if ord==Inf, opts.quad = 'm'; end
  
  Ns=50:50:450; err = zeros(size(Ns)); tic;
  for i=1:numel(Ns);
    N = Ns(i);
    s = scale(segment.smoothstar(N, 0.3, 5), 0.5);  % same pentafoil as in pc2d
    S = layerpot.S(k, s, [], opts);
    D = layerpot.D(k, s, [], opts);
    if k==0
      f = real(s.x).*imag(s.x);           % bdry vals quadrupole field u=xy
      g = imag(s.x).*real(s.nx) + real(s.x).*imag(s.nx);% grad u =(y,x)
    else
      c = exp(1i*pi/3);                  % bdry vals u = plane wave, direc c
      f = exp(1i*k*real(conj(c)*s.x));
      g = 1i*k*real(conj(c)*s.nx).*f;  % grad u = ik<c.normal>*value
    end
    tau = -f; sig = g;   % Greens Rep Thm choice for densities
    err(i) = norm((D - eye(N)/2)*tau + S*sig - f)*sqrt(2*pi/N); % L_2 on bdry
    disp(sprintf('N=%d\tL_2 err = %.16g', N, err(i)));
  end
  toc
  loglog(Ns, err, ['+-' co(o)]);
  hold on; plot(Ns, 1e-2*(Ns/1e2).^-ord, ['--' co(o)]); drawnow
end % -----------------------------------------
axis tight; legend(num2cellstr(kron(ords, [1 1])));
xlabel('N'); ylabel('L_2 err')
title([sprintf('k=%g inside Greens Rep convergence at orders: ',k) ...
       sprintf('%d ', ords)]);
if verb>1, print -depsc2 test_lpquad_GRF_helmholtz.eps
end
if verb, figure; imagesc(real(S)); colorbar; end



% ============================== crude test Alpert vs Kress
N = 50; k=10;
s = segment.smoothstar(N, 0.3, 3);
opts.quad = 'k';
SK = layerpot.S(k, s, [], opts);
opts.quad = 'a';
S = layerpot.S(k, s, [], opts);
norm(S-SK)/norm(S)
