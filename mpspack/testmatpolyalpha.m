% build up and test alpha matrix polynomial for polyeig
% barnett 8/20/08

clear classes
test = 'sl';    % first char: '1' one beta value, 's' beta sweep, 'd' a,b sweep
                %   'p' for (a,b) sweep using alpha-polynomial
                % 2nd char: QP basis: 'j' for J-exp, 'l' for LPs
inv = 1;        % 0 compute whole BZ, 1 use inversion symm
sys = 'd';      % system: 'e' for empty unit cell, 'd' for single dirichlet blob
verb = 1;         % plot type (eg eigs for test 's')
k = 12; trgnei = 0; srcnei = 0;
Mu = 12;
uc = qpunitcell(1, 0.5+1i, k, Mu);   % M=20 for uc quadr ('g' default type)
ns = 100;                       % max # sing vals to keep
Ms = 30; s = scale(segment.smoothstar(Ms, 0.3, 3), 0.2);
e = domain([], [], s, -1);
b = e.addlayerpotbasis(s, 'd', k);  % replace 'd' by [1i*k 1] for BWLP
b = b{1}; % get basis handle in b

square = 0;
if test(2)=='j'
  o.rescale_rad = 1.0; uc.addregfbbasis(0, 20, k, o); % N=5 for k=2
elseif test(2)=='l'
  uc.addqpuclayerpots; square = 1;
end

if test(1)=='1'   % ..................................... basic tests
  if sys=='e', be = exp(1i*pi*0.58); else be = 1; end % (b ang = 0.58, Bindel)
  uc.setbloch(exp(1i*pi*0.82638948), be);
  o.poly = 1;
  tic; Q = uc.evalbasesdiscrep(o); fprintf('Q from cold: %.2g s\n', toc)
  tic; Q = uc.evalbasesdiscrep(o); fprintf('Q reusing data: %.2g s\n', toc)
  if sys=='d'
    o.nei = trgnei; tic; [B dB] = uc.evalbasestargetcopies(s, o); toc
    o.nei = srcnei; o.dom = e; tic; [C dC] = b.evalunitcelldiscrep(uc, o); toc
    o.sourcenei = srcnei; o.targetnei = trgnei;
    tic; [A dA] = b.evalsourcecopiestargetcopies(s, uc, o); toc
    M = [2*A 2*B; C Q];  % 2's correspond to halving defn of mismatch
  else, M = Q; end
  if verb>1, figure; showmatpoly(M, 'M'); end
  if ~square         % econ size QR gets random onb in P...
    [P dummy] = qr(rand(size(M,1), size(M,2)),0);
    Mo = M; M = zeros(size(Mo,2),size(Mo,2),size(Mo,3));
    for i=1:5, M(:,:,i) = P'*Mo(:,:,i); end  % project M to make sq (Bindel idea)
    end
  fprintf('polyeig...\n');
  tic; E = polyeig(squeeze(M(:,:,1)), squeeze(M(:,:,2)), squeeze(M(:,:,3)), squeeze(M(:,:,4)), squeeze(M(:,:,5))); fprintf('polyeig %.2g s\n', toc)
  if verb, figure; plot(E, '+'); axis equal; axis([-2 2 -2 2]);title('eigs');
    hold on; circle(0,0,1,'k-'); end
  alphas = E(find(abs(abs(E)-1)<1e-3));        % the alpha vals found by polyeig
  aangs = pi*(-1:0.02:1);     % test these eigs against alpha angle sweep
  o.poly = 0; ss = NaN*zeros(ns, numel(aangs)); tic;
  for i=1:numel(aangs)
    uc.setbloch(exp(1i*aangs(i)), be);
    Q = uc.evalbasesdiscrep(o);
    if sys=='d', o.nei = trgnei; o.data = dB;
      B = uc.evalbasestargetcopies(s, o); % dB is data
      o.nei = srcnei; o.dom = e; o.data = dC;
      C = b.evalunitcelldiscrep(uc, o);
      o.sourcenei = srcnei; o.targetnei = trgnei; o.data = dA;
      A = b.evalsourcecopiestargetcopies(s, uc, o);
      M = [2*A 2*B; C Q];  % 2's correspond to halving defn of mismatch
    else, M = Q; end
    sing = svd(M);
    n = min(ns,numel(sing)); ss(1:n,i) = sing(end:-1:end-n+1);
  end
  fprintf('alpha sweep %.2g s\n', toc)
  if verb, figure; plot(aangs, ss, '.-'); hold on;
    if ~isempty(alphas), vline(angle(alphas)); end
    xlabel('ang \alpha'); ylabel('sing val of M'); v=axis; v(4)=1; axis(v); end
  
  if verb>1,        % do some relevant timing tests....
    tic; N=320; eig(rand(N), rand(N)); toc % about half the polyeig time at N=80
    tic; for i=1:100, svd(M); end; toc
  end
    
  
elseif test(1)=='s'   % ...................................... beta sweep of BZ
  % verb choices: 1 BZ plot, 2 alpha eig movie, 3 log alpha eig movie.
  o.poly = 1; bangs = pi*(0:0.05:1);  % beta angles, new:(0.42:0.001:0.48)@k=10
  if verb==1, bf=figure; uc.setbloch(1, 1); uc.showbrillouin; hold on; R = uc.recip/(2*pi); end
  if verb>1, ef=figure; end
  uc.setbloch(1, 1);
  if sys=='d'
    tic; Q = uc.evalbasesdiscrep(o); toc
    o.nei = trgnei; tic; [B dB] = uc.evalbasestargetcopies(s, o); toc
    o.nei = srcnei; o.dom = e; %o = rmfield(o, 'data');
    tic; [C dC] = b.evalunitcelldiscrep(uc, o); toc
    o.sourcenei = srcnei; o.targetnei = trgnei;
    tic; [A dA] = b.evalsourcecopiestargetcopies(s, uc, o); toc
    M = [2*A 2*B; C Q];  % 2's correspond to halving defn of mismatch
  end
  tic;
  for i=1:numel(bangs)
    uc.setbloch(1, exp(1i*bangs(i)));
    Q = uc.evalbasesdiscrep(o);
    if sys=='d', o.nei = trgnei; o.data =dB; B = uc.evalbasestargetcopies(s, o);
      o.nei = srcnei; o.dom = e; o.data = dC; C = b.evalunitcelldiscrep(uc, o);
      o.sourcenei = srcnei; o.targetnei = trgnei; o.data = dA;
      A = b.evalsourcecopiestargetcopies(s, uc, o);
      M = [2*A 2*B; C Q];  % 2's correspond to halving defn of mismatch
    else, M = Q; end
    if ~square         % econ size QR gets random onb in P...
      [P dummy] = qr(rand(size(M,1), size(M,2)),0);
      Mo = M; M = zeros(size(Mo,2),size(Mo,2),size(Mo,3));
      for j=1:5, M(:,:,j) = P'*Mo(:,:,j); end
    end
    E = polyeig(squeeze(M(:,:,1)), squeeze(M(:,:,2)), squeeze(M(:,:,3)), squeeze(M(:,:,4)));       % notice dropped the a^2 term since zero for nei=0
    alphas = E(find(abs(abs(E)-1)<1e-3));   % the alpha vals found by polyeig
    if verb==1
      if inv
        c = R*[angle([alphas; conj(alphas)])'; ...
               bangs(i)*kron([1 -1], ones(1,numel(alphas)))];
      else c = R*[angle(alphas)'; bangs(i)*ones(1,numel(alphas))]; end
      figure(bf); plot(c(1,:), c(2,:), '.'); axis equal; end
    if verb==2, figure(ef);plot(E(find(abs(log10(abs(E)))<3)), '+');axis equal;
      axis([-2 2 -2 2]); circle(0,0,1,'k-'); title('\alpha eigs'); end
    if verb==3, figure(ef); plot(log(E(find(abs(log10(abs(E)))<3))), '+');
      axis([-pi pi -3 3]); vline(0,'k-'); title('log \alpha eigs'); end
    if verb>1,text(1,1,sprintf('\\beta ang/\\pi=%.4g',bangs(i)/pi));drawnow; end
  end
  toc
  if verb==1, xlabel('Bloch k_x'); ylabel('Bloch k_y');
    set(gcf, 'name', 'Brillouin zone'); text(pi,pi,sprintf('k=%g',k));
    if sys=='e', circle(0,0,k,'k-'); end               % k-circle
  end
  
  
elseif test(1)=='d'   % ....................... alpha,beta double sweep, k=const
  o.poly = 0;
  if sys=='d' % fill dA dB dC for use in loop...
    o.nei = trgnei; tic; [B dB] = uc.evalbasestargetcopies(s, o); toc
    o.nei = srcnei; o.dom = e; tic; [C dC] = b.evalunitcelldiscrep(uc, o); toc
    o.sourcenei = srcnei; o.targetnei = trgnei;
    tic; [A dA] = b.evalsourcecopiestargetcopies(s, uc, o); toc
  end
  aangs = pi*(-1:0.05:1);    % alpha angles
  bangs = pi*(0:0.05:1);     % beta angles
  ss = NaN*zeros(ns, numel(aangs), numel(bangs));
  profile clear; profile on;
  for i=1:numel(aangs)
    fprintf('i=%d\n', i)
    for j=1:numel(bangs)
      uc.setbloch(exp(1i*aangs(i)), exp(1i*bangs(j)));
      Q = uc.evalbasesdiscrep(o);
      if sys=='d', o.nei = trgnei; o.data=dB; B=uc.evalbasestargetcopies(s, o);
        o.nei = srcnei; o.dom=e; o.data=dC; C = b.evalunitcelldiscrep(uc, o);
        o.sourcenei = srcnei; o.targetnei = trgnei; o.data = dA;
        A = b.evalsourcecopiestargetcopies(s, uc, o);
        M = [2*A 2*B; C Q];  % 2's correspond to halving defn of mismatch
      else, M = Q; end
      sing = svd(M);
      n = min(ns,numel(sing)); ss(1:n,i,j) = sing(end:-1:end-n+1);
    end
  end
  profile off; profile viewer
  if verb, figure; uc.setbloch(1, 1);                               % BZ plot
    uc.showbrillouin; hold on; R = uc.recip/(2*pi);
    [aa bb] = meshgrid(aangs, bangs); ab = R*[aa(:)'; bb(:)'];
    aa = reshape(ab(1,:), size(aa)); bb = reshape(ab(2,:), size(aa));
    surf(aa, bb, 0*aa, log10(squeeze(ss(1,:,:))'));
    colorbar; shading interp; axis equal; set(gcf, 'name', 'Brillouin zone');
    set(gca, 'ydir', 'normal'); xlabel('Bloch k_x'); ylabel('Bloch k_y');
    if sys=='e', circle(0,0,k,'k-'); end              % k-circle
  end
  if verb>1, figure; imagesc(aangs, bangs, log10(squeeze(ss(1,:,:)))); colorbar;
    set(gca, 'ydir', 'normal'); xlabel('ang \alpha'); ylabel('ang \beta');
    figure; plot(aangs, squeeze(ss(:,11,:)), '+-'); 
  end
  
  
elseif test(1)=='p'   % .... use alpha-poly for alpha,beta double sweep, k=const
  o.poly = 1;
  if sys=='d' % fill dA dB dC for use in loop...
    o.nei = trgnei; tic; [B dB] = uc.evalbasestargetcopies(s, o); toc
    o.nei = srcnei; o.dom = e; tic; [C dC] = b.evalunitcelldiscrep(uc, o); toc
    o.sourcenei = srcnei; o.targetnei = trgnei;
    tic; [A dA] = b.evalsourcecopiestargetcopies(s, uc, o); toc
  end
  aangs = pi*(-1:0.05:1);     % alpha angles
  bangs = pi*(-1:0.05:1);     % beta angles
  if inv, bangs = pi*(0:0.05:1); end % bottom half is inversion symm, don't eval
  ss = NaN*zeros(numel(aangs), numel(bangs), ns);
  tic; %profile clear; profile on;
  for j=1:numel(bangs)
    uc.setbloch(1, exp(1i*bangs(j)));
    %fprintf('j=%d\tbeta ang=%.3f pi\n',j,bangs(j)/pi)    
    Q = uc.evalbasesdiscrep(o);              % set up alpha-poly at fixed beta
    if sys=='d', o.nei = trgnei; o.data=dB; B=uc.evalbasestargetcopies(s, o);
      o.nei = srcnei; o.dom=e; o.data=dC; C = b.evalunitcelldiscrep(uc, o);
      o.sourcenei = srcnei; o.targetnei = trgnei; o.data = dA;
      A = b.evalsourcecopiestargetcopies(s, uc, o);
      M = [2*A 2*B; C Q];  % 2's correspond to halving defn of mismatch
    else, M = Q; end
    v = ones(size(M,2),1);            % for use by lin solve
                                      %v2 = reshape(1:numel(v), size(v));
    for i=1:numel(aangs)
      al = exp(1i*aangs(i));
      Ma = M(:,:,1)*al^-2 + M(:,:,2)*al^-1 + M(:,:,3) + M(:,:,4)*al + M(:,:,5)*al^2;           % evaluate the matrix poly 
      %sing = svd(Ma); n=min(ns,numel(sing)); ss(i,j,1:n)=sing(end:-1:end-n+1);
      ss(i,j,1) = 1/norm(Ma \ v);
      %ss(i,j,1) = 1/(norm(Ma \ v) + norm(Ma \ v2));  % attempt at more robust?
    end
  end
  toc %profile off; profile viewer
  if inv
    ss = [ss(end:-1:1,end:-1:2,:) ss]; bangs=[-bangs(end:-1:2) bangs]; % invsymm
  end
  if verb, figure; uc.setbloch(1, 1);                               % BZ plot
    uc.showbrillouin; hold on; R = uc.recip/(2*pi);
    [aa bb] = meshgrid(aangs, bangs); ab = R*[aa(:)'; bb(:)'];
    aa = reshape(ab(1,:), size(aa)); bb = reshape(ab(2,:), size(aa));
    surf(aa, bb, 0*aa, log10(squeeze(ss(:,:,1))'));
    colorbar; shading interp; axis equal; set(gcf, 'name', 'Brillouin zone');
    set(gca, 'ydir', 'normal'); xlabel('Bloch k_x'); ylabel('Bloch k_y');
    v = axis;
    if sys=='e', for i=-2:2, for j=-2:2, x = i*uc.r1+j*uc.r2;
          circle(real(x),imag(x),k,'k-'); end, end, end              % k-circles
    axis(v);
  end
  if verb>1, figure; imagesc(aangs, bangs, log10(squeeze(ss(:,:,1)))); colorbar;
    set(gca, 'ydir', 'normal'); xlabel('ang \alpha'); ylabel('ang \beta');
    figure; plot(aangs, squeeze(ss(:,11,:)), '+-'); 
  end
end
