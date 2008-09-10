% Bloch mode problem without the 'problem' class, for single Dirichlet segment.
% barnett 8/22/08, cleaned up and changed for simpler copylists, 9/9/08

clear all classes
sys = 'e';             % 'e' empty UC, 's' scatterer (Dirichlet closed seg)
bas = 'q';             % QP basis type ('j' is only old test; 'q' is real deal)
% note test cases specific to k=10, eg sys='s' reaches 1e-14 for Mu=30,Ms=70
k = 10; Mu = 30; Ms = 46;                % Mu quadr pts per UC wall seg
verb = 0; bwcan = 1;             % verbosity, & whether to cancel Bloch wave
uc = qpunitcell(1, 0.5+1i, k, Mu);       % decide the UC shape (1, 0.5+1i)
uc.buffer = 1;       % choose QP scheme. 0: 1xUC sticking out (nei=0,1)
                     % 1: 3xUC sticking out (nei=2), 2: 3x scaled UC (nei=1,2)
o.nei = 2;           % 0,1,etc: # segment src neighbor copies (1x1, 3x3, etc)
if uc.buffer==2, ucQsiz = 3;                    % NB won't work for alpha-poly
  ucQ = qpunitcell(3*uc.e1, 3*uc.e2, k, Mu);    % discrep uses 3x size unit cell
else ucQ = uc; ucQsiz = 1; end                  % discrep uses original uc
if bas=='j', ucQ.addregfbbasis(0,15,k); uc.setupbasisdofs;
elseif bas=='q', ucQ.addqpuclayerpots; end
s = scale(segment.smoothstar(Ms, 0.3, 3), 0.2);  % Ms quadr pts on inclusion
e = domain([], [], s, -1); o.dom = e;
b = e.addlayerpotbasis(s, [1i*k 1], k); b = b{1}; % BWLP basis; get basis handle
if verb>1, figure; uc.plot; hold on; ucQ.plot; s.plot; title('uc & ucQ'); end

if 1 % no poly...
  if sys=='s', uc.setbloch(exp(1i*pi*0.82638945966117), 1); % a mode to 12 digs
  else uc.setbloch(exp(1i*pi*0.818199331592963),1i); end    % a mode to 15 digs
  ucQ.setbloch(uc.a^ucQsiz, uc.b^ucQsiz);   % maybe cubed discrepancy phases
  fprintf('Filling times (Q,B,C,A) from cold:\n');
  tic; Q = ucQ.evalbasesdiscrep; toc % poly data storage is in ucQ automatically
  if sys~='e', tic; [B dB] = ucQ.evalbaseswithdata(s); toc % (ucQ>uc poly bad)
  tic; [C dC] = b.evalunitcelldiscrep(uc, o); toc % heeds ucbuf=2
  tic; [A dA] = b.evalunitcellcopies(s, uc, o); toc % plain square source block
  M = [2*A 2*B; C Q]; if verb>1,figure;imagesc(real(M));colorbar;title('M');end
  else M = Q; A = []; end
  [U S V] = svd(M); sig = diag(S); svec = V(:,end); % min right sing vec=coeffs
  fprintf('min sing val of M = %.15g   (should be v small)\n', sig(end))
  if verb, if uc.buffer~=1, g = -1:0.01:1; else g = -3:0.03:3; end % grid mode
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

if 0 % no poly... using stored values, check code speed profile
  if sys=='s', uc.setbloch(exp(1i*pi*-0.82638945966), 1);  % also mode
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
  
if 0 % poly... from stored data, then do a PEP matrix solve on it...
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