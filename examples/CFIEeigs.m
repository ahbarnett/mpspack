% expts on eigs of A_k = 1/2 + D - i.k.S for CFIE, w/ Euan Spence
% Barnett 3/6/19

clear all classes; verb = 1;          % if verb>0, generates EPS figures
if 0   % smooth
  tref = segment.radialfunc([], {@(q) 1 + 0.3*cos(3*q), @(q) -0.9*sin(3*q),...
                      @(q) -2.7*cos(3*q)});
  zs = 0.2+0.3i;
  Ns = 10:10:100;
  Ns = 200:200:600;
elseif 0      % corner cavity
  o.kressq=8;
  %tref = segment.polyseglist(50,[-1-1i, 1-1i, 1+1i, -1+1i],'pc',o);
  %ss = [-1-1i, 1-1i, 1-.5i, -.5i, .5i, 1+.5i, 1+1i, -1+1i];
  a = 0.2; ss = [-1-1i, 1-1i, 1-.3i, -a-.5i-1i*a, -a+.5i+1i*a, 1+.3i, 1+1i, -1+1i];  % cavity
  tref = segment.polyseglist(50,ss,'pc',o);
  zs = -0.5+.2i;
Ns = 10:10:100;
else                        % smooth resonant cup
  tref = larrycup(100);                    % *** change to new seg type, move
                                % larrycup out of exmaples/
  zs = 1.0+0.1i;
%  Ns = 40:40:500;   % sufficient for k=40;   280 enough for k=20 
Ns = 200:200:600;   % ok for k=80
end
d = domain([], [], tref, -1);
%d.k = 40; % 6*pi; %20;
%d.k = 18.550;   % res for cup
%d.k = 10.68;
%d.k = 14.62;
%d.k = fzero(@(x) besselj(0,x),15)/0.8;  % J_0 mode k for 0.8 inner rad
%d.k=19.56;  % 1.2-scaled cup res
%d.k=39.1735;  % 1.2-scaled cup res, bowtie odd
%d.k=39.89;  % 1.2-scaled cup res, bowtie even
%d.k=40.43;  % 1.2-scaled cup messy res
%d.k=41.136;  % 1.2-scaled cup: 2nd-hermite laser cavity mode
%d.k=41.87; %42.39093;
d.k = 80.22;

figure; d.plot; drawnow
d.addlayerpot(tref, 'D');                    % adds DLP to tref segment
f = @(z) besselh(0,d.k * abs(z-zs)); % known exterior field
for i=1:numel(tref)
  tref(i).setbc(1, 'D', [], @(t) f(tref(i).Z(t))); % its Dirichlet data
end
p = bvp(d);
for i=1:numel(Ns)
  p.updateN(Ns(i)); p.solvecoeffs; N(i) = p.N;
  e(i) = abs(f(2) - p.pointsolution(pointset(2)));  % single pt soln test
end
cond(p.A)
%figure; p.showsolution; drawnow
p.solvecoeffs(struct('meth','iter','eps',1e-10));

%figure; semilogy(N, e, '+-'); xlabel('N'); ylabel('error in u(2)');
if verb,       % generate f.lp a,b
  g=figure; set(gca,'fontsize', 20); semilogy(N, e, '+-'); axis tight;
  xlabel('N'); ylabel('abs error in u(2)');
  set(gcf,'paperposition', [.25 .25 6 8]);
  print -depsc2 ../doc/figs/lpconv.eps
  figure; plot(eig(diag(1./p.sqrtwei)*p.A), '+'); set(gca,'fontsize', 20);
  axis([-.1 1.1 -.6 .6]); axis equal;
  hold on; t=0:0.01:2*pi; plot(0.5 + 0.5*exp(1i*t), 'r-');
  xlabel('Re[\lambda(1/2+D)]'); ylabel('Im[\lambda(1/2+D)]');
  set(gcf,'paperposition', [.25 .25 6 8]);
%  print -depsc2 ../doc/figs/lpeig.eps
end

% Demo BWLP combined-field... (note I changed -ikS sign on Timo's suggestion)
p.bas{1}.a = [-1i*d.k 1];           % sneaky way to change SLP,DLP coeffs
p.fillbcmatrix;
cond(p.A)
p.solvecoeffs(struct('meth','iter','eps',1e-10));
if verb,       % regenerate f.lp a, and generate f.lp c
  for i=1:numel(Ns)
    p.updateN(Ns(i)); p.solvecoeffs; N(i) = p.N;
    e(i) = abs(f(2) - p.pointsolution(pointset(2)));
  end
  figure(g); hold on; semilogy(N, e, 'go--'); axis tight;
  print -depsc2 ../doc/figs/lpconv.eps
  figure; plot(eig(diag(1./p.sqrtwei)*p.A), '+'); set(gca,'fontsize', 20);
  hold on; t=0:0.01:2*pi; plot(0.5 + 0.5*exp(1i*t), 'r-');
  axis equal tight; axis(1.05*axis);
  xlabel('Re[\lambda(1/2+D-ikS)]'); ylabel('Im[\lambda(1/2+D-ikS)]');
  set(gcf,'paperposition', [.25 .25 6 8]);
%  print -depsc2 ../doc/figs/lpeig_bwlp.eps
end
figure; p.showsolution; drawnow
figure; plot(eig(diag(1./p.sqrtwei)*p.A), '+'); axis equal
hold on; t=0:0.01:2*pi; plot(0.5 + 0.5*exp(1i*t), 'r-'); hold off
    
[V D] = eig(diag(1./p.sqrtwei)*p.A);
j=find(abs(diag(D))==min(abs(diag(D))));
fprintf('min eig(A_k) = %.3g + %.3g i\n',real(D(j,j)),imag(D(j,j)))
tau =V(:,j);
p.co = tau;
figure; p.showsolution(struct('dx',0.02)); drawnow; title('mode from min eigval of A_k');

stop

%%%%%%%%%%%%%%%%%%%%%%%  k anim, for CFIE (fixed eta=k)
figure; %p.updateN(70);
ks = 80:0.01:81; %18:0.02:21; %ks = (5.5:0.01:6.5)*pi;
p.bas{1}.a = [-1i*d.k 1];           % sneaky way to change SLP,DLP coeffs
for i=1:numel(ks)
  d.k = ks(i);
  p.fillbcmatrix;
  plot(eig(diag(1./p.sqrtwei)*p.A), '+'); axis equal
  axis([-.5 3.5 -2 2]);
    hold on; t=0:0.01:2*pi; plot(0.5 + 0.5*exp(1i*t), 'r-'); hold off
    %title(sprintf('k = %.3g pi \n',ks(i)/pi))
    title(sprintf('k = %.3f\n',ks(i)))
  drawnow;
end





