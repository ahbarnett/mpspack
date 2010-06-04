% Example codes from MPSpack tutorial, also generates EPS figures for this doc
% SECTION 6: LAYER POTENTIALS

clear all classes; verb = 1;          % if verb>0, generates EPS figures
tref = segment.radialfunc([], {@(q) 1 + 0.3*cos(3*q), @(q) -0.9*sin(3*q),...
                   @(q) -2.7*cos(3*q)});
d = domain([], [], tref, -1); d.k = 10;
d.addlayerpot(tref, 'D');                    % adds DLP to tref segment
f = @(z) besselh(0,d.k * abs(z-0.3-0.2i)); % known exterior field
tref.setbc(1, 'D', [], @(t) f(tref.Z(t))); % its Dirichlet data
p = bvp(d);
Ns = 5:5:80; for i=1:numel(Ns)
  p.updateN(Ns(i)); p.solvecoeffs; N(i) = p.N;
  e(i) = abs(f(2) - p.pointsolution(pointset(2)));
end
cond(p.A)
figure; semilogy(N, e, '+-'); xlabel('N'); ylabel('error in u(2)');
if verb,       % generate f.lp a,b
  g=figure; set(gca,'fontsize', 20); semilogy(N, e, '+-'); axis tight;
  xlabel('N'); ylabel('abs error in u(2)');
  set(gcf,'paperposition', [.25 .25 6 8]);
  print -depsc2 ../doc/lpconv.eps
  figure; plot(eig(diag(1./p.sqrtwei)*p.A), '+'); set(gca,'fontsize', 20);
  axis([-.1 1.1 -.6 .6]); axis equal;
  hold on; t=0:0.01:2*pi; plot(0.5 + 0.5*exp(1i*t), 'r-');
  xlabel('Re[\lambda(1/2+D)]'); ylabel('Im[\lambda(1/2+D)]');
  set(gcf,'paperposition', [.25 .25 6 8]);
  print -depsc2 ../doc/lpeig.eps
end

% Demo BWLP combined-field... (note I changed -ikS sign on Timo's suggestion)
p.bas{1}.a = [-1i*d.k 1];           % sneaky way to change SLP,DLP coeffs
p.fillbcmatrix;
cond(p.A)
if verb,       % regenerate f.lp a, and generate f.lp c
  Ns = 5:5:80; for i=1:numel(Ns)
    p.updateN(Ns(i)); p.solvecoeffs; N(i) = p.N;
    e(i) = abs(f(2) - p.pointsolution(pointset(2)));
  end
  figure(g); hold on; semilogy(N, e, 'go--'); axis tight;
  print -depsc2 ../doc/lpconv.eps
  figure; plot(eig(diag(1./p.sqrtwei)*p.A), '+'); set(gca,'fontsize', 20);
  hold on; t=0:0.01:2*pi; plot(0.5 + 0.5*exp(1i*t), 'r-');
  axis equal tight; axis(1.05*axis);
  xlabel('Re[\lambda(1/2+D-ikS)]'); ylabel('Im[\lambda(1/2+D-ikS)]');
  set(gcf,'paperposition', [.25 .25 6 8]);
  print -depsc2 ../doc/lpeig_bwlp.eps
end
