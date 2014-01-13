% test Kapur-Rokhlin, Kress, quadrature on single & double layer potentials
% barnett 2/5/08 adapted to MPSpack 7/31/08.
% This is lower-level test routine than testlayerpot
clear classes
verb = 0;            % verbosity

if 0 % =============== test SLP/DLP vals, and in Laplace case const DLP Greens
k = 0;                      % wavenumber (>=0 for 's' or 'd', =0 for 'c')
ords = [2 6 10 Inf 2 6 10];  % Kapur orders (Inf gives Kress, -ve Alpert)
schemes = 'kkkmaaa'; co = 'bgrkcym';         % quad schemes & plot colors
test = 'c';                  % 'd'/'s' for DLP/SLP vals, 'c' for const DLP k=0

if test=='c', figure; end
for o = 1:numel(ords) % -------------------------- loop over orders
  opts.ord = ords(o); opts.quad = schemes(o);
  
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
    hold on; plot(Ns, 1e-2*(Ns/1e2).^-opts.ord, ['--' co(o)]); drawnow % slope
  end
end % -----------------------------------------
if test=='c', axis tight; legend(num2cellstr(kron(ords, [1 1])));
  xlabel('N'); ylabel('L_2 err');
  title(['DLP k=0 tau\equiv1 convergence at orders: ' sprintf('%d ', ords)]);
  if verb>1, print -depsc2 test_dlpquad_laplace.eps
  end
end
if verb, figure; imagesc(real(D)); colorbar; figure; show_bdry(bd); end
end

if 1 % ==================== Test Greens Thm for combined SLP/DLP, k>=0
k = 10;                      % wavenumber (may be 0 or >0)
ords = [2 6 10 Inf 2 6 10 16 pi]; % orders (Inf for Kress)
schemes = 'kkkmaaaaa'; co = 'bgrkcymbg';         % quad schemes & plot colors

figure;
for o = 1:numel(ords) % -------------------------- loop over orders
  opts.ord = ords(o); opts.quad = schemes(o);
  
  Ns=[50 100 200 300 400 500 600]; err = zeros(size(Ns)); tic;
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
  loglog(Ns, err, ['+-' co(o)]); axis([min(Ns) max(Ns) 1e-16 1]);
  hold on; %plot(Ns, 1e-2*(Ns/1e2).^-opts.ord, ['--' co(o)]);
  drawnow
end % -----------------------------------------
%axis tight; legend(num2cellstr(kron(ords, [1 1])),'location','southwest');
axis tight; legend(num2cellstr(ords),'location','southwest');
xlabel('N'); ylabel('L_2 err')
title([sprintf('k=%g inside Greens Rep convergence at orders: ',k) ...
       sprintf('%d ', ords)]);
if verb>1, print -depsc2 lpquad_GRF_helmholtz_Alpert+2ninterp.eps
end
if verb, figure; imagesc(real(S)); colorbar; end
end


if 0 % ============================ warm-up: crude test Alpert vs Kress for SLP
clear all classed; N = 50; k=10;
s = segment.smoothstar(N, 0.3, 3);
opts.quad = 'm';  % Kress
SK = layerpot.S(k, s, [], opts);
figure; imagesc(real(SK)); colorbar; caxis(1e-1*[-1 1]);
opts.quad = 'a'; opts.ord = 6; % Alpert
profile clear; profile on;
S = layerpot.S(k, s, [], opts);
profile off; profile viewer;
figure; imagesc(real(S)); colorbar; caxis(1e-1*[-1 1]); % looks similar!
sigma = cos(2*pi*(1:N)'/N);
f = S*sigma; fk = SK*sigma;       % compare their effect on smooth density
figure; plot((1:N)/N, real([f fk]), '+-');  % good!
end

if 0 % ================================== timing tests for Alpertizing things
  clear all classes; s = segment.smoothstar(450, 0.3, 3); k = 10;
  opts.quad = 'a'; opts.ord = 6; % Alpert
  profile clear; profile on;
  tic; S = layerpot.S(k, s, [], opts); toc
  profile off; profile viewer;

  tic; n=1e4; for i=1:n, s.Z(rand(1)); end; toc  % 600 us per eval!
  tic; a = rand(1,n); s.Z(a); toc            % 0.6 us per eval!
  tic; for i=1:n, besselh(0,rand(1)); end; toc  % 120 us per eval!
  tic; a = rand(1,n); besselh(0,a); toc            % 2 us per eval!  
end