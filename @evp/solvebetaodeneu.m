function kh = solvebetaodeneu(kstar, beo, A, B, C, E)
% SOLVEBETAODENEU - helper for higher-order khat Neumann DtN scaling method
%
% kh = SOLVEBETAODE(kstar, beo, A, B, C, E) solves the ODE
%   db/dk = -kA + (B - Cb - Eb^2)/k,    with IC b(kstar) = beo,
%   returning the wavenumber value kh at which b(kh) = 0.
%
% Experimental!!! based on solvebetaode, 1/15/14

% x = kstar - k, is indep var of ODE (more accurate than using k):
F = @(x,b) +(kstar-x)*A - (B - C*b - E.*b.*b)./(kstar-x);  % use k -> kstar-x
op = odeset('abstol',1e-14,'events',@ev); %'initialstep',0.5*kstar*abs(beo));
[x,b,xev,bev,iev] = ode45(F, [0 1.0], real(beo), op); % 1.0>eps
kh = kstar - xev; % k for the event happening (stopping-point)
%figure; plot(x,b,'+-'); % debug
if isempty(kh), warning(sprintf('solvebetaodeneu failed to find intersection event! beo=%g x(end)=%g b(end)=%g\n', beo, x(end), b(end)));
  kh = kstar + beo/kstar; % back to standard 1st-order approx
end
  
function [v,ist,dir] = ev(x,b)  % event func for ode45 stopping
v = b; ist = 1; dir = +1; % if b were complex, need v=real(b)
