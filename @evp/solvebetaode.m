function kh = solvebetaode(kstar, beo, A, B, C, Ap, Bp, Cp)
% SOLVEBETAODE - helper routine for higher-order khat in NtD scaling method
%
% kh = SOLVEBETAODE(kstar, beo, A, B, C) solves the ODE
%   db/dk = (1 - Bb + (k^2.C - A)b^2)/k,    with IC b(kstar) = beo,
%   returning the wavenumber value kh at which b(kh) = 0.
%
%   This is used to give a better khat prediction.
%   kstar = starting k, beo = starting beta, A=||x_t f_t||^2, B = Re<f,mf>,
%   C = ||x_t f||^2. kh = predicted \hat{k}
%
% kh = SOLVEBETAODE(kstar, beo, A, B, C, Ap, Bp, Cp) uses a linear approximation
%   instead of frozen values. Ap,Bp,Cp are derivatives A,B,C taken wrt k.
%   This seems to offer little improvement.
%
% Notes:
% Andrew's fixed ODE, 6/12/11. Alex's attempt to use lin approx for f, 6/13/11

% x = k-kstar is indep var of ODE (more accurate than using k):
if nargin==5   % ODE for b(x), fixed Andrew's...
  F = @(x,b) (1 - B*b + ((kstar+x).^2*C-A).*b.*b)./(kstar+x);
elseif nargin==8
  F = @(x,b) (1 - (B + Bp*x).*b + ((kstar+x).^2*(C+Cp*x)-(A+Ap*x)).*b.*b)./(kstar+x); % Alex
end
op = odeset('abstol',1e-14,'events',@ev); %'initialstep',0.5*kstar*abs(beo));
[x,b,xev,bev,iev] = ode45(F, [0 1.0], real(beo), op); % 1.0>eps
kh = xev + kstar; % k for the event happening (stopping-point)
%figure; plot(x,b,'+-'); % debug
if isempty(kh), warning(sprintf('solvebetaode failed to find intersection event! beo=%g x(end)=%g b(end)=%g\n', beo, x(end), b(end)));
  kh = kstar/(1+beo); % back to standard 1st-order approx
end
  
function [v,ist,dir] = ev(x,b)  % event func for ode45 stopping
v = b; ist = 1; dir = +1; % if b were complex, need v=real(b)
