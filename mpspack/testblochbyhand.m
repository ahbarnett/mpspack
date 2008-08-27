% bloch mode problem without the 'problem' class, for single Dirichlet segment
% barnett 8/22/08

clear all classes
k = 10; Mu = 20; Ms = 50; verb = 1;
srcnei = 0; trgnei = 0; % size of source, target neighborhoods (1x1 or 3x3)
uc = qpunitcell(1, 0.5+1i, k, Mu);
uc.setbloch(exp(1i*pi*0.82638948), 1);  % pi*0.82638948 is close to mode
uc.addqpuclayerpots;
s = scale(segment.smoothstar(Ms, 0.3, 3), 0.2);
e = domain([], [], s, -1);
b = e.addlayerpotbasis(s, 'd', k); b = b{1}; % get basis handle
if verb>1, figure; uc.plot; hold on; s.plot; end

if 0 % nail down B problem ! Bpoly output shoudn't dep on alpha in uc!
  o.poly = 1; uc.setbloch(1, 1);
  tic; [B dB] = uc.evalbasestargetcopies(s, o); toc % dB is data
  uc.setbloch(1i, 1); %o.data = dB;
  tic; B2 = uc.evalbasestargetcopies(s, o); toc % dB is data
  figure; showmatpoly(B2 - B, 'B2-B'); o.poly = 0;
end

if 1 % no poly...
  uc.setbloch(exp(1i*pi*0.82638948), 1);  % pi*0.82638948 is close to mode
  tic; Q = uc.evalbasesdiscrep; toc % poly data storage is in uc automatically
  o.nei = trgnei; tic; [B dB] = uc.evalbasestargetcopies(s, o); toc % dB is data
  o.nei = srcnei; o.dom = e; tic; [C dC] = b.evalunitcelldiscrep(uc, o); toc
  o.sourcenei = srcnei; o.targetnei = trgnei;
  tic; [A dA] = b.evalsourcecopiestargetcopies(s, uc, o); toc
  M = [2*A 2*B; C Q]; %figure; imagesc(real(M)); colorbar; title('M');
  min(svd(M)) % should be about <1e-8
end

if 1 % poly...
o.poly = 1; o.dom = e; tic; Q = uc.evalbasesdiscrep(o); toc
o.nei = trgnei; o.data = dB;
tic; [B dB] = uc.evalbasestargetcopies(s, o); toc % dB is data
save B B
o.nei = srcnei; o.dom = e; o.data = dC;
tic; [C dC] = b.evalunitcelldiscrep(uc, o); toc
o.sourcenei = srcnei; o.targetnei = trgnei; o.data = dA;
tic; [A dA] = b.evalsourcecopiestargetcopies(s, uc, o); toc

M = [2*A 2*B; C Q];  % 2's correspond to halving defn of mismatch (to get Id)
if verb, figure; showmatpoly(M,'M'); end
if verb & exist('Mcold'), figure; showmatpoly(M-Mcold,'M-Mcold'); end
fprintf('polyeig...\n');
tic; E = polyeig(squeeze(M(:,:,1)), squeeze(M(:,:,2)), squeeze(M(:,:,3)), squeeze(M(:,:,4)), squeeze(M(:,:,5))); fprintf('polyeig %.2g s\n', toc)
if verb, figure; plot(E, '+'); axis equal; axis([-2 2 -2 2]);title('eigs');
  hold on; circle(0,0,1,'k-'); end
alphas = E(find(abs(abs(E)-1)<1e-3));        % the alpha vals found by polyeig
angle(alphas)/pi
abs(alphas)
end
