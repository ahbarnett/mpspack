% test ability to compute discrepancy for generic & QPLP basis functions
% barnett 9/7/08

clear classes
k = 9; o.ucbuf = 0; o.nei = 0; bas = 'j';  % basis type 'j' or 'q' (QPLPs)
verb = 1;         % plot verbosity
inv = 1;          % 0 compute whole BZ, 1 use inversion symm
Mu = 12;
uc = qpunitcell(1, 0.5+1i, k, Mu);   % M=20 for uc quadr ('g' default type)
if bas=='j'
  o.rescale_rad = 1.0; uc.addregfbbasis(0, 10, k, o); % N=5 for k=2, N=10 k=9
elseif bas=='l'
  uc.addqpuclayerpots;
end
uc.setupbasisdofs;
o.poly = 0; ns = 100;
aangs = pi*(-1:0.05:1);    % alpha angles
bangs = pi*(-1:0.05:1);     % beta angles
if inv, bangs = pi*(0:0.05:1); end % bottom half is inversion symm, don't eval
ss = NaN*zeros(numel(aangs), numel(bangs), ns);
%profile clear; profile on;
 for i=1:numel(aangs)
   fprintf('i=%d\n', i)
   for j=1:numel(bangs)
     uc.setbloch(exp(1i*aangs(i)), exp(1i*bangs(j)));
     %[Q o.data] = uc.bas{1}.evalunitcelldiscrep(uc, o);
     Q = uc.evalbasesdiscrep(o);
     sing = svd(Q); n = min(ns,numel(sing)); ss(i,j,1:n) = sing(end:-1:end-n+1);
   end
  end
%profile off; profile viewer
if inv
  ss = [ss(end:-1:1,end:-1:2,:) ss]; bangs=[-bangs(end:-1:2) bangs]; % invsymm
end
if verb, figure; uc.setbloch(1, 1);                               % BZ plot
  uc.showbrillouin; hold on; R = uc.recip/(2*pi);
  [aa bb] = meshgrid(aangs, bangs); ab = R*[aa(:)'; bb(:)'];
  aa = reshape(ab(1,:), size(aa)); bb = reshape(ab(2,:), size(aa));
  surf(aa, bb, 0*aa, log10(squeeze(ss(:,:,1))'));
  colorbar; shading interp; axis equal; set(gcf, 'name', 'Brillouin zone');
  set(gca, 'ydir', 'normal'); xlabel('Bloch k_x'); ylabel('Bloch k_y'); v=axis;
  for i=-2:2, for j=-2:2, x = i*uc.r1+j*uc.r2;
      circle(real(x),imag(x),k,'k-'); end, end; axis(v);          % k-circles
end
