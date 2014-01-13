% test Alpert log-quadratures on periodic log singularity's Fourier coeffs
% barnett 1/29/11, added Alpert's osc*log quadr 12/16/13
clear
ords = [2,3,4,5,6,8,10,12,14,16,pi,Inf];   % order "pi" is 2013 new quadr
P = 2*pi;     % period
m = 10;   % Kress periodized log * m^th Fourier mode (incr m for more challenge)
f = @(t) log(4*sin(t/2).^2) .* cos(2*pi*m*t/P)/P;
Iex = 0; if m~=0, Iex = -1/abs(m); end   % exact answer (Kress LIE book)
Ns = 10:2:70;    % #s of periodic nodes to try (even for Kress to work)
IN = nan(numel(Ns), numel(ords));

for j=1:numel(ords), ord = ords(j); fprintf('\n ord = %d:\n', ord);
  for i=1:numel(Ns), N = Ns(i); h = P/N;
    if isinf(ord)      % check Kress for kicks (f func is log sing times sth)
      t = h*(0:N-1); w = quadr.kress_Rjn(N/2);  % h is already in w
      IN(i,j) = sum(w .* cos(2*pi*m*t/P)/P);  % explicitly removed f's log sing
    else               % Alpert
      [tex,wex,nskip] = quadr.QuadLogExtraPtNodes(ord); % extra nodes
      if N>2*nskip
        t = [tex; (nskip:(N-nskip))'; N-tex(end:-1:1)] * h; % nodes
        w = [wex; ones(N-2*nskip+1, 1); wex(end:-1:1)] * h; % weights
        IN(i,j) = sum(f(t).*w);
      end
    end
    fprintf('N=%d: err = %g\n', N, abs(IN(i,j)-Iex))
  end
end
figure; loglog(Ns, abs(IN-Iex), '+-'); legnum(ords); hold on;
plot(Ns, repmat(Ns.'/m, [1 numel(ords)]).^-repmat(ords, [numel(Ns) 1]), '--');
xlabel('N'); ylabel('error'); title('Alpert log endpoint correction errors');
axis([min(Ns) max(Ns) 1e-16 1]);
