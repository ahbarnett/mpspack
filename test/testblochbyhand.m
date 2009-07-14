% Bloch mode problem without the 'problem' class, for single segment inclusion.
% Barnett 8/22/08. cleaned up and changed for simpler copylists, 9/9/08
% Jfilter 1/20/09, cleaned up for diel transmission inclusion 1/23/09

clear all classes
sys = 't';  % 'e' empty UC, 's' Dirichlet scatt, 'n' Neumann, 't' transm
n=0.3; k = 14.9; Mu = 35; Ms = 80;
%n=3; k = 7.8; Mu = 20; Ms = 50;       % Mu quadr pts per UC wall seg, eg 40 for Jfilt
%n = 2; k = 5; Mu = 20; Ms = 40;
verb = 1; opts.verb = verb; opts.bwcan = 0; % verbosity & if cancel Bloch wave
uc = qpunitcell(1, 0.5+1i, k, Mu);  % UC lattice vectors (1, 0.5+1i)
uc.buffer = 0; o.nei = 1; % choose QP scheme. 0: 1xUC sticking out (nei=0,1)
bo=[]; bo.Jfilter.M = 16; bo.Jfilter.rescale_rad = 0.5;% B: J-filter QPLP eval
uc.addqpuclayerpots;
s = scale(segment.smoothstar(Ms, 0.3, 3), 0.42); % 0.42, Ms quadr pts on inclusion
%s.translate(0.1-.35i);      % test ew invariance under transl (>=0.23 hits UC)
de = domain([], [], s, -1);   % infinite domain exterior to inclusion
if sys=='s'
  b = de.addlayerpotbasis(s, [1i*k 1], k); b = b{1}; % BWLP pot get basis handle
  uc.setbloch(exp(1i*pi*0.82638945966117), 1); % a mode to 12 digs
elseif sys=='n'
  b = de.addlayerpotbasis(s, 's', k); b = b{1}; % SLP (no BWLP since cot absent)
  uc.setbloch(exp(1i*pi*0.4675075815), 1); % a mode
elseif sys=='t'
  di = domain(s, 1);          % interior domain of inclusion
  de.addlayerpotbasis(s, 's', k); de.addlayerpotbasis(s, 'd', k); % NB order!
  di.addlayerpotbasis(s, 's', n*k); di.addlayerpotbasis(s, 'd', n*k); % "
  b = [de.bas{:} di.bas{:}];          % two ext and two int LP bases = 4 total
  uc.setbloch(exp(1i*pi*-0.125084764525647),exp(1i*pi*0.6)); %for: k5 s0.2 n2
  uc.setbloch(exp(1i*pi*0.0557817719182246),exp(1i*pi*0.7)); % k14.9 s0.42 n0.3
%  uc.setbloch(exp(1i*pi*0.0557817719182246),exp(1i*pi*0)); % k14.9 s0.42 n0.3
elseif sys=='e'
  b = [];                     % dummy arg to fillblochmodematrix
  uc.setbloch(exp(1i*pi*0.818199331592963),1i); % empty mode to 15 digs (k=14.9)
end
if verb>2, figure; uc.plot; hold on; s.plot; title('UC & inclusion'); end

if 1 % no poly...
  [M data] = fillblochmodematrix(sys, b, uc, s, o, bo, opts); opts.data = data; 
  if verb>2,figure;imagesc(real(M));colorbar;title('M');end
  [U S V] = svd(M); sig = diag(S); svec = V(:,end); % min right sing vec=coeffs
  fprintf('min sing val of M = %.15g   (should be v small)\n', sig(end))
  if verb>1, opts.skew=1; if ~isempty(bo), opts.Jfilter=bo.Jfilter; end
    showblochmode(sys, b, uc, s, svec, o, opts); end
  if sys=='e' & verb>1          % test empty mode against a plane wave
    bwave = exp(-1i*real(conj(uc.kbloch-4*pi).*(xx+1i*yy))); % Bloch * plane wv
    c = u(101,101).*bwave(101,101);  % value at center of grid, test u=const?
    imagesc(g, g, log10(abs((u+v).*bwave/c-1))); set(gca, 'ydir','normal');
    colorbar; uc.plot; axis equal;
    title('relative error in empty Bloch mode, k=10, Mu=30');
  end
  if verb>2    % movie to show complex values by cycling phase before Re[] ...
    f=figure; N=300; for j=1:N, ph=exp(1i*2*pi*j/N); figure(f);
      imagesc(g, g, real(ph*(u+v).*bwave)); set(gca, 'ydir','normal');
      caxis(.02*[-1 1]); axis equal; drawnow; end
  end
end

if 0   % ------------------------------- sweep alpha & beta over BZ
  da = 0.03; inv = 1;       % # singvals keep, & whether to use BZ inv symm
  aangs = angle(uc.a); %pi*(-1:da:1);     % alpha angles (make eg 0 for single slice)
  bangs = pi*(-1:da:1);     % beta angles
  if inv, bangs = pi*(0:da:1); end % bottom half is inversion symm, don't eval
  ss = NaN*zeros(numel(aangs), numel(bangs), min(size(M))); opts.verb = 0;
  t = timetic; tic(t); rhsvec = rand(size(M,1),1) - 0.5;
  for i=1:numel(aangs)
    fprintf('i=%d/%d\n', i, numel(aangs))
    for j=1:numel(bangs)
      uc.setbloch(exp(1i*aangs(i)), exp(1i*bangs(j)));
      M = fillblochmodematrix(sys, b, uc, s, o, bo, opts);
      sing = svd(M); ss(i,j,:) = sing(end:-1:1);
      ss(i,j,1) = 1./norm(M\rhsvec);
    end
  end
  fprintf('Brillouin zone sweep (%d SVDs) done in %g s\n', numel(aangs)*numel(bangs), toc(t))
  if inv
    ss = [ss(end:-1:1,end:-1:2,:) ss]; bangs=[-bangs(end:-1:2) bangs]; % invsymm
  end
  if numel(aangs)==1
    figure; plot(bangs, squeeze(ss(1,:,1:10)), '-'); xlabel('\beta ang');
    axis([0 3 0 1e-1]);
  else
  figure; uc.setbloch(1, 1); uc.showbrillouin; hold on; R = uc.recip/(2*pi);
  [aa bb] = meshgrid(aangs, bangs); ab = R*[aa(:)'; bb(:)'];
  aa = reshape(ab(1,:), size(aa)); bb = reshape(ab(2,:), size(aa));
  surf(aa, bb, 0*aa, log10(squeeze(ss(:,:,1))'));
  colorbar; shading interp; axis equal;
  set(gcf, 'name', sprintf('k=%g: Brillouin zone', k));
  set(gca, 'ydir', 'normal'); xlabel('Bloch k_x'); ylabel('Bloch k_y');
  v = axis;
  if sys=='e', for i=-2:2, for j=-2:2, x = i*uc.r1+j*uc.r2;
        circle(real(x),imag(x),k,'k-'); end, end, end              % k-circles
        axis(v);
  end
  %save diel_s0.2_n=3_k7.8_slice.mat
end


if 0 % no poly... using stored values, check code speed profile
  if sys=='s', uc.setbloch(exp(1i*pi*-0.82638945966), 1);  % also mode
  else uc.setbloch(exp(1i*pi*-0.618810864777384),1i); end % also empty mode
  fprintf('Filling times (Q,B,C,A) reusing stored data, at a new k_block:\n');
  %profile clear; profile on; for i=1:100
  M = fillblochmodematrix(sys, b, uc, s, o, bo, opts);
  %end; profile off; profile viewer
  fprintf('min sing val of M = %.15g     (should also be small)\n', min(svd(M)))
end
  
if 1 % poly... from stored data, then do a PEP matrix solve on it...
  o.poly = 3; bo.poly = o.poly; % eg cubic for 1xUC scheme
  fprintf('Filling (Q,B,C,A) with poly=%d, stored data:\n', o.poly);
  [M opts.data] = fillblochmodematrix(sys, b, uc, s, o, bo, opts);
  if verb>1, figure; showmatpoly(M,'M'); end
  fprintf('solving PEP via polyeig (linearized dense GEP order N=%d)...\n',...
          size(M,1)*o.poly);
  Mcells = num2cell(M, [1 2]); % break into cell array along dimension 3 only
  tic; E = polyeig(Mcells{:}); % break cells into comma-separated list
    fprintf('polyeig %.2g s\n', toc)
  if verb, figure; plot(E, '+'); axis equal; axis([-2 2 -2 2]);title('eigs');
    hold on; circle(0,0,1,'k-'); end
  alphas = E(find(abs(abs(E)-1)<1e-3)); % keep PEP evals lying near unit circle
  if numel(alphas)>1, if sys=='s'
    fprintf('eigval angs/pi: %.15g %.15g   (close to +-0.82638946)\n', ...
            angle(alphas(1))/pi, angle(alphas(2))/pi)
    fprintf('|eigvals|-1:    %g %g   (should be small)\n', ...
            abs(alphas(1))-1, abs(alphas(2))-1)
  else
    disp('eigval angles/pi :'); angle(alphas)/pi
    disp('|eigvals|-1 :'); abs(alphas)-1
    end
  o.poly = 0; bo.poly = 0; j=1;  % choose which poly-found eigenmode to plot...
  uc.a = exp(1i*angle(alphas(j)));
  M = fillblochmodematrix(sys, b, uc, s, o, bo, opts);
  [U S V] = svd(M); sig = diag(S); svec = V(:,end); % min right sing vec=coeffs
  fprintf('min sing val of M = %.15g   (should be v small)\n', sig(end))
  opts.uctrim = 0; opts.gy = -.6:0.01:.6; opts.gx= -.8:0.01:.8;
  if ~isempty(bo), opts.Jfilter=bo.Jfilter; end 
if verb, [u gx gy]=showblochmode(sys, b, uc, s, svec, o, opts); end
  %save diel_s0.42_n0.3_k9_b0.6pi_j1_mode.mat
  else, disp('no eigvals found'); end
end

if 0 % debug ucbuf=1
  ucb = 1;
  figure; uc.bas{ucb}.showQdatasegs(uc.Qpolydata{ucb}); uc.plot; axis equal
end

if 0   % speed tests
  profile clear; profile on
  for i=1:100
    uc = qpunitcell(1, 0.5+1i, k, 2); uc.addqpuclayerpots; % NB only Mu=2!
    Q = uc.evalbasesdiscrep;
  end
  profile off; profile viewer
end
  
if 0   % ------------------------------- volume sweep omega, alpha & beta
  ks =  0.1:0.03:10;
  nsv = min(10, size(M,1));   % how many singvals to keep (< dim M)
  na = 20; inv = 1; % # samples per pi unit, whether to use BZ inv symm
  aangs = pi*((1:2*na)/na-1);   % alpha angles
  bangs = pi*((1:2*na)/na-1);    % beta angles
  ss = NaN*zeros(numel(ks), numel(aangs), numel(bangs), nsv);
  opts.verb = 0;
  for w=1:numel(ks)
    k = ks(w); fprintf('k=%g:\n', k); t1 = timetic; tic(t1);
    uc.k = k; uc.Qpolydata=[]; for i=1:4,uc.bas{i}.k=k; end % change k (messy!)
    if sys=='s' | sys=='n', de.bas{1}.k = k;
    elseif sys=='t'
      de.bas{1}.k = k; de.bas{2}.k = k; di.bas{1}.k = n*k; di.bas{2}.k = n*k; 
    end
    [M opts.data] = fillblochmodematrix(sys, b, uc, s, o, bo); % change k !
    for i=1:numel(aangs)
      fprintf('\ti=%d/%d\n', i, numel(aangs))
      for j=1:numel(bangs)
        if j<=na | j==2*na
          uc.setbloch(exp(1i*aangs(i)), exp(1i*bangs(j)));   % calculate
          M = fillblochmodematrix(sys, b, uc, s, o, bo, opts);
          sing = svd(M);  %[U S V] = svd(M); sing = diag(S); % twice as slow
          ss(w,i,j,:) = sing(end:-1:end-nsv+1);
        end
      end
    end
    fprintf('\tk done in %g s\n', toc(t1))
  end
  if inv,       % use inversion k-space symm (t-rev invariance)
    ss(:, 2*na, na+1:2*na-1, :) = ss(:, 2*na, na-1:-1:1, :);  % inv symm X-M
    ss(:,1:2*na-1,na+1:2*na-1,:) = ss(:,2*na-1:-1:1,na-1:-1:1,:);
  end
  nam = 'diel_s0.2_n=3_bandvolume'; 
  %save(sprintf('%s.mat',nam));
  % concatenate then plot a journey around Brioullin zone.
  slic = [squeeze(ss(:,na+1:end,na,1)) squeeze(ss(:,end,na+1:end,1))];
  for i=1:na-1, slic = [slic squeeze(ss(:,end-i,end-i,1))]; end
  figure; sp = abs(uc.r1+uc.r2)/2/pi;  % relative speed along brill diag 
  set(gca,'position', [.1 .1 .4 .8], 'fontsize', 14);
  imagesc([0:na, na+sp*(1:na) na*(1+sp)+(1:na)], ks, log10(slic));
  caxis([-4 -1]); colormap(gray(256)); set(gca, 'ydir','normal', 'xtick', []);
  ylabel('\omega'); vline(na, 'k-');vline(2*na, 'k-'); vline((2+sp)*na, 'k-');
  vline(0,'k-'); text(0,-0.3,'\Gamma');text(na,-0.3,'X');text(2*na,-0.3,'M');
  text((2+sp)*na,-0.3,'\Gamma');
  print('-depsc2', sprintf('%s_slicetour.eps',nam));
end

if 0  % .............................. higher sing vals at const a,b
  figure; plot(ks, squeeze(ss(:,21,21,:)), '-');
end

if 0 % .......... make movie of above data
figure; uc.setbloch(1, 1); R = uc.recip/(2*pi);
[aa bb] = meshgrid(aangs, bangs); ab = R*[aa(:)'; bb(:)'];
aa = reshape(ab(1,:), size(aa)); bb = reshape(ab(2,:), size(aa));
set(gcf, 'position', [300 300 300 240]); clear Mov
for j=1:numel(ks),
  uc.showbrillouin; hold on;
  surf(aa, bb, 0*aa, log10(squeeze(ss(j,:,:,1))'));
  %colorbar;
  shading interp; axis equal;
  xlabel('Bloch k_x'); ylabel('Bloch k_y');
  set(gca, 'ydir','normal'); text(1,3,sprintf('k=%.2f',ks(j)));
  axis equal; axis tight;
  caxis([-6 -1]); colormap(jet(256)); %caxis([-4 0]); %colormap(hot(256));
  drawnow;
  Mov(j) = getframe(gcf); %pause(.1);
  clf;
end
movie2avi(Mov,sprintf('%s.avi',nam));   % writes uncompressed AVI
system(sprintf('mencoder %s.avi -o %s.mpg -ovc lavc', nam, nam)); % make mpeg
% write a mini gif for webpage... (note -zoom needed for software zoom via -xy)
system(sprintf('mplayer %s.avi -zoom -xy 0.5 -vo gif89a:fps=15.0:output=%s.gif', nam, nam));
end
