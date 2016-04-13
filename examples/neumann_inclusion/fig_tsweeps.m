% tension sweep subfigs for B-Hassell-Tacy paper. Barnett 12/8/15
clear

krng=1;
if krng==1 % low k:
  r = load('rf_neu_ref40k41_N450_p_bdry'); kjex = r.p.kj; jdot = 12; % blob ind
  load rfn_newtsweep_40k41
  kslvalid = [min(kjex) max(kjex)];      % interval of k to check slopes in
  sc = 1;
else
  kjex = 405.0032695182281;   jdot = 1;
  load('rfn_newtsweep_405k405.1.mat');
  kslvalid = [405.0003 405.0057];
  sc = 1e3;    % hack so ticks never have exponent which screws EPS bbox!
end  
tmins = min(ts,[],1);
figure; plot(kswp.^2/sc,tmins,'-'); %plot(kswp.^2,ts,'b.'); % all curves
if krng==2
  set(gca,'xticklabel',[164040:20:164100],'xtick',[164040:20:164100]/sc);
end
xlabel('$E$','interpreter','latex');
ylabel('$\tilde{t}_{h,min}(E)$','interpreter','latex');
axis([min(kswp)^2/sc max(kswp)^2/sc 0 3.2])
% vline(kjex.^2);  % not needed
hold on; plot(kjex(jdot)^2/sc, 0, 'r.','markersize',20)

if krng==1, text(min(kswp)^2+3,2.7,'(a)');
set(gcf,'paperposition',[0 0 5.1 1.5]);
  print -depsc2 tsweep_med.eps
else,       text((min(kswp)^2+3)/sc,2.7,'(b)');
  set(gcf,'paperposition',[0 0 5 1.5]);
  print -depsc2 tsweep_high.eps
end
  
% measure E slopes of tension
Ejex = kjex.^2; ne = numel(kjex);
nearEs = min(abs(repmat(kswp.^2,[ne 1])-repmat(Ejex,[1 numel(kswp)])),[],1);
ii = kswp>kslvalid(1) & kswp<kslvalid(2) & nearEs>0.01;  % valid for slopes
sl = tmins(ii)./nearEs(ii); % slopes for not-too-close cases
fprintf('E tension slope range = [%.3g, %.3g]\n',min(sl),max(sl))
%[0.646, 0.676]
Cest = 1/min(sl)  % estimate for const in our theorem: 1.6


