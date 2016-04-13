% accessibility figure 1 for neumann paper, 10/28/15
k = 100;
h = 1/k;
chi1 = @(t) (1+erf(5*(t-1.5)))/2;  % effectively smooth cutoff
%chi2 = @(t) sqrt(1-chi1(t).^2);
%figure; t = 0:0.01:3; plot(t,[chi1(t);chi2(t)], '-'); title('\chi_1, \chi_2')
chi2 = @(t) 1 - chi1(t);
posqrt = @(s) sqrt(s).*(s>=0);
Gh = @(s) sqrt(s).*chi1(s/h^(2/3)) + h^(1/3)*chi2(s/h^(2/3));
smax = 5.0*h^(2/3); s = smax*(-.5:0.01:1.5);
figure; subplot(2,1,1); plot(s,posqrt(s), 'k--', 'linewidth', 2); hold on;
plot(s,Gh(s),'k-', 'linewidth', 2); axis tight; xlabel('\sigma');
legend('(\sigma)^{1/2}_+','G_h(\sigma)','location','northeast');
text(-0.4*smax,.5,'(c)');
x = -1.3:0.01:1.7;
subplot(2,1,2); plot(x,posqrt(1-x.^2), 'k--', 'linewidth', 2); hold on;
plot(x,Gh(1-x.^2), 'k-', 'linewidth', 2);
xlabel('\xi'''); axis tight;
legend('(1-\xi''^2)^{1/2}_+','G_h(1-\xi''^2)');
text(-1.2,.8,'(d)');
set(gcf,'paperposition',[0 0 5 4]);
print -depsc2 cutoff.eps
