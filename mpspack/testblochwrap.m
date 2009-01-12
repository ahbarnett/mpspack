% Bloch mode problem without the 'problem' class, for single Dirichlet segment.
% barnett 8/22/08, cleaned up and changed for simpler copylists, 9/9/08
% wrapping of 1xUC scheme (R wall only; including poly), 1/11/09.

clear all classes
sys = 's';             % 'e' empty UC, 's' scatterer (Dirichlet closed seg)
bas = 'q';             % QP basis type ('j' is only test; 'q' is real deal)
% note test cases specific to k=10, eg sys='s' reaches 1e-14 for Mu=30,Ms=70
k = 10; Mu = 20; Ms = 40;           % Mu quadr pts per UC wall seg (18,80)
verb = 0; bwcan = 1;             % verbosity, & whether to cancel Bloch wave
uc = qpunitcell(1, 0.5+1i, k, Mu);       % decide the UC shape (1, 0.5+1i)
uc.buffer = 0;       % choose QP scheme. 0: 1xUC sticking out (nei=0,1)
                     % 1: 3xUC sticking out (nei=2), 2: 3x scaled UC (nei=1,2)
o.nei = 1;           % 0,1,etc: # segment src neighbor copies (1x1, 3x3, etc)
o.close = 0.25;   % dist to go adaptive for close-evaluation in B filling
if uc.buffer==2, ucQsiz = 3;                    % NB won't work for alpha-poly
  ucQ = qpunitcell(3*uc.e1, 3*uc.e2, k, Mu);    % discrep uses 3x size unit cell
else ucQ = uc; ucQsiz = 1; end                  % discrep uses original uc
if bas=='j', ucQ.addregfbbasis(0,38,k); ucQ.setupbasisdofs;
elseif bas=='q', ucQ.addqpuclayerpots; end
s = scale(segment.smoothstar(Ms, 0.3, 3), 0.2);  % Ms quadr pts on inclusion
%s = scale(segment.smoothstar(Ms, 0.3, 3), 0.42); s.translate(0.1);
s.translate(0.3);        % test ew invariance under transl (>0.229 hits UC R)
wrap = 1;                % 0: none, 1: wrap A & B correctly (R wall only)
if wrap==1 & (o.nei~=1 | uc.buffer~=0), disp('o.nei must be 1 to wrap!'); end
e = domain([], [], s, -1); o.dom = e;
b = e.addlayerpotbasis(s, [1i*k 1], k); b = b{1}; % BWLP pot; get basis handle
%b = e.addlayerpotbasis(s, 's', k); b = b{1}; % BWLP pot; get basis handle
if verb>1, figure; uc.plot; hold on; ucQ.plot; s.plot; title('uc & ucQ'); end

if 1 % no poly...
  if sys=='s', uc.setbloch(exp(1i*pi*0.82638945966117), 1); % a mode to 12 digs
  else uc.setbloch(exp(1i*pi*0.818199331592963),1i); end % empty mode to 15 digs
%uc.setbloch(exp(1i*pi*(0.818199331592963-2/3)),1i);  % spurious 3xempty eigval
  ucQ.setbloch(uc.a^ucQsiz, uc.b^ucQsiz);   % maybe cubed discrepancy phases
  fprintf('Filling times (Q,B,A,C) from cold:\n');
  tic; Q = ucQ.evalbasesdiscrep; toc % poly data storage is in ucQ automatically
  if sys~='e'
    if wrap==1, oldsx = s.x; [s.x ix jx] = uc.fold(s.x); end
    tic; [B dB] = ucQ.evalbaseswithdata(s,o); toc % (ucQ>uc poly bad)
    if wrap==1, B = diag(uc.a.^ix.*uc.b.^jx)*B; s.x=oldsx; % rephase B rows
      for i=1:numel(dB), dB{i} = uc.datawrapR(dB{i}, find(ix>0)); end, end
    o.close=0;                                     % don't use close for A, C
    tic; [A dA] = b.evalunitcellcopies(s, uc, o); toc % square source block
    if wrap==1, ii=find(ix>0);  % wrapped rows, spec to G2 penetrating R UC wall
      cc = find(dA.copylist.apow==-1); dA.B(ii,:,cc) = 0; %kill (-1,:) copy rows
      o.copylist.t = -2*uc.e1 - (-1:1)*uc.e2; % new (2,:) col of copies for G2
      o.copylist.apow = [2 2 2]; o.copylist.bpow = [-1 0 1]; % (above trans -ve)
      o.copylist.remph = [1 1 1]; p2 = pointset(s.x(ii), s.nx(ii)); % G2 pts
      [A2 dA2] = b.evalunitcellcopies(p2, uc, o); % new copies (2,:)
      ns=size(dA.B,3)+(1:3); dA.B(:,:,ns) = 0; dA.B(ii,:,ns) = dA2.B;
      dA.copylist.t(ns) = o.copylist.t; dA.copylist.apow(ns) = o.copylist.apow;
      dA.copylist.bpow(ns) = o.copylist.bpow; dA.copylist.remph(ns) = 1;
      o = rmfield(o, 'copylist'); clear A2 dA2;
      o.data = dA; A = b.evalunitcellcopies(s, uc, o); % resum w/ new list
      o = rmfield(o, 'data');   % now dA has correct data for this freq k
      end
    tic; [C dC] = b.evalunitcelldiscrep(uc, o); toc % heeds ucbuf=2
    M = [2*A 2*B; C Q];
    if verb>1,figure;imagesc(real(M));colorbar;title('M');end
  else M = Q; A = []; end
  [U S V] = svd(M); sig = diag(S); svec = V(:,end); % min right sing vec=coeffs
  fprintf('min sing val of M = %.15g   (should be v small)\n', sig(end))
  if verb, if uc.buffer==0, g = -1:0.01:1; else g = -2:0.02:2; end % grid mode
    [xx yy]=meshgrid(g,g); p=pointset(xx(:)+1i*yy(:)); disp('eval QP field...')
    tic;B = ucQ.evalbases(p);toc,u = reshape(B*svec(size(A,2)+1:end), size(xx));
    if sys~='e', disp('eval obst field...'); tic;
      A = b.evalunitcellcopies(p, uc, o); toc
      v = reshape(A*svec(1:size(A,2)), size(xx)); % u is QP field, v obst field
    else v=0; end
    bwave = 1; if bwcan,
      bwave = exp(-1i*real(conj(uc.kbloch).*(xx+1i*yy))); end % conj bloch wave
    figure; imagesc(g, g, real((u+v).*bwave)); set(gca, 'ydir','normal');
    uc.plot; xlabel x; ylabel y; title('QP + obst field w/ Bloch wave cancel');
    axis equal; caxis(.02*[-1 1]);
  end
  if sys=='e' & verb>1          % test empty mode against a plane wave
    bwave = exp(-1i*real(conj(uc.kbloch-4*pi).*(xx+1i*yy))); % Bloch * plane wv
    c = u(101,101).*bwave(101,101);  % value at center of grid, test u=const?
    imagesc(g, g, log10(abs((u+v).*bwave/c-1))); set(gca, 'ydir','normal');
    colorbar; uc.plot; axis equal;
    title('relative error in empty Bloch mode, k=10, ucbuf=2, Mu=30');
    print -depsc2 empty_k10_ucbuf2.eps
  end
  if verb>1    % movie to show complex values by cycling phase before Re[] ...
    f=figure; N=300; for j=1:N, ph=exp(1i*2*pi*j/N); figure(f);
      imagesc(g, g, real(ph*(u+v).*bwave)); set(gca, 'ydir','normal');
      caxis(.02*[-1 1]); axis equal; drawnow; end
  end
end

if 1            % sweep alpha & beta over BZ
  da = 0.1; inv = 1; % # singvals keep, & whether to use BZ inv symm
  aangs = pi*(-1:da:1);    % alpha angles  (make 0 for single slice)
  bangs = pi*(-1:da:1);     % beta angles
  %aangs = 2.5:0.01:2.7; bangs = -0.3:0.01:-0.1; weak area for restr-u-nrm trick
  %bangs = 2.7:1e-3:2.8;
  %bangs = 2.726155:1e-8:2.726157;     % beta angles 2.72615612 is mode for a=0
  if inv, bangs = pi*(0:da:1); end % bottom half is inversion symm, don't eval
  ss = NaN*zeros(numel(aangs), numel(bangs), min(size(M))); sss = ss;
  v = zeros(size(M,2),1); v(1:size(A,2)) = 1;  % restricted RHS vector
  for i=1:numel(aangs)
    fprintf('i=%d\n', i)
    for j=1:numel(bangs)
      uc.setbloch(exp(1i*aangs(i)), exp(1i*bangs(j)));
      ucQ.setbloch(uc.a^ucQsiz, uc.b^ucQsiz);  % maybe cubed discrepancy phases
      Q = ucQ.evalbasesdiscrep;
      if sys~='e', o.data = dB; B = ucQ.evalbaseswithdata(s, o);
        o.data = dC; C = b.evalunitcelldiscrep(uc, o);
        o.data = dA; A = b.evalunitcellcopies(s, uc, o);
        M = [2*A 2*B; C Q]; else M = Q; end
      if uc.buffer>0
        [U S V] = svd(M); sing = diag(S); sss(i,j,:) = sing(end:-1:1); % incr ord
        ss(i,j,:)=sing(end:-1:1)'./(uc.spuriousdistance+sqrt(sum(abs(U(1:size(A,2),end:-1:1)).^2,1)));
      %ss(i,j,end:-1:1) = sum(abs(U(1:size(A,2),:)).^2,1); % restric u nrm^2
      else, sing = svd(M); ss(i,j,:) = sing(end:-1:1); end       % default
      %ss(i,j,:) = ss(i,j,:) / uc.spuriousdistance;      % kill known spurious
      %ews = eig(M); ss(i,j,:) = sort(abs(ews), 'ascend');  % try eigs
      %ss(i,j,1) = 1/norm(M\v);       % use restric-rhs lin solve (gmres?)
    end
  end
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
  colorbar; shading interp; axis equal; set(gcf, 'name', 'Brillouin zone');
  set(gca, 'ydir', 'normal'); xlabel('Bloch k_x'); ylabel('Bloch k_y');
  v = axis;
  if sys=='e', for i=-2:2, for j=-2:2, x = i*uc.r1+j*uc.r2;
        circle(real(x),imag(x),k,'k-'); end, end, end              % k-circles
 % if sys=='s', for i=-2:1/3:2, for j=-3:1/3:3, x = i*uc.r1+j*uc.r2;
 %       if floor(i)-i~=0 | floor(j)-j~=0, circle(real(x),imag(x),k,'k-'); end
 %     end, end, end              % spurious k-circles
    axis(v);
  end
end

if 1 % no poly... using stored values, check code speed profile
  if sys=='s', uc.setbloch(exp(1i*pi*-0.82638945966117), 1);  % also mode
  else uc.setbloch(exp(1i*pi*-0.618810864777384),1i); end % also empty mode
  fprintf('Filling times (Q,B,C,A) reusing stored data, at a new k_block:\n');
       %profile clear; profile on; for i=1:100
  tic; Q = uc.evalbasesdiscrep; toc % poly data storage is in uc automatically
  if sys~='e', tic; o.data = dB; B = uc.evalbaseswithdata(s, o); toc
  tic; o.data = dC; C = b.evalunitcelldiscrep(uc, o); toc
  tic; o.data = dA; A = b.evalunitcellcopies(s, uc, o); toc
       %end; profile off; profile viewer
  M = [2*A 2*B; C Q];
  else M = Q; end
  fprintf('min sing val of M = %.15g     (should also be small)\n', min(svd(M)))
end
  
if 0 % poly... from stored data, then do a PEP matrix solve on it... wrap ok too
  o.poly = 4 + 3*uc.buffer;  % quartic, or septic for full scheme
  fprintf('Filling (Q,B,C,A) with poly=%d, stored data:\n', o.poly);
  tic; Q = uc.evalbasesdiscrep(o); toc % poly data store is in uc automatically
  if sys~='e',tic; o.data = dB; B = uc.evalbaseswithdata(s, o); toc % dB is data
  tic; o.data = dC; C = b.evalunitcelldiscrep(uc, o); toc
  tic; o.data = dA; A = b.evalunitcellcopies(s, uc, o); toc
  M = [2*A 2*B; C Q];  % 2's correspond to halving defn of mismatch (to get Id)
  else M=Q; end
  if verb, figure; showmatpoly(M,'M'); end
  fprintf('solving PEP via polyeig (linearized dense GEP order N=%d)...\n',...
          size(M,1)*o.poly);
  Mcells = num2cell(M, [1 2]); % break into cell array along dimension 3 only
  tic; E = polyeig(Mcells{:}); % break cells into comma-separated list
    fprintf('polyeig %.2g s\n', toc)
  if verb, figure; plot(E, '+'); axis equal; axis([-2 2 -2 2]);title('eigs');
    hold on; circle(0,0,1,'k-'); end
  alphas = E(find(abs(abs(E)-1)<1e-3)); % keep PEP evals lying near unit circle
  if sys=='s'
    fprintf('eigval angs/pi: %.15g %.15g   (close to +-0.82638946)\n', ...
            angle(alphas(1))/pi, angle(alphas(2))/pi)
    fprintf('|eigvals|-1:    %g %g   (should be small)\n', ...
            abs(alphas(1))-1, abs(alphas(2))-1)
  else
    disp('eigval angles/pi :'); angle(alphas)/pi
    disp('|eigvals|-1 :'); abs(alphas)-1
  end
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