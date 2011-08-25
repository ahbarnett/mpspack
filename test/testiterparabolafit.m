% test the EVP local minimization routine against MATLAB's fminbnd
% Barnett 8/23/11

clear; x = [-1 0 1];  % used throughout

f = @(x) x^2;          % the vanilla minimum
[xm fm] = fminbnd(f, x(1), x(3), optimset('display','iter'))  % matlab
o.xtol = 1e-12; o.verb = 1; for i=1:3, y(i) = f(x(i)); end  % setup input
[xm fm] = evp.iterparabolafit(f, x, y, o)               % mine

f = @(x) (exp(x)-1)^2;   % skew minimum
[xm fm] = fminbnd(f, x(1), x(3), optimset('display','iter'))
o.xtol = 1e-12; o.verb = 1; for i=1:3, y(i) = f(x(i)); end
[xm fm i] = evp.iterparabolafit(f, x, y, o)

f = @(x) (exp(x))^2;   % minimum beyond left of interval
[xm fm] = fminbnd(f, x(1), x(3), optimset('display','iter'))
o.xtol = 1e-12; o.verb = 1; for i=1:3, y(i) = f(x(i)); end
[xm fm i] = evp.iterparabolafit(f, x, y, o)

f = @(x) -x + 1e-8*x^2;   % minimum to right of interval, weak curve
[xm fm] = fminbnd(f, x(1), x(3), optimset('display','iter'))
o.xtol = 1e-12; o.verb = 1; for i=1:3, y(i) = f(x(i)); end
[xm fm] = evp.iterparabolafit(f, x, y, o)

f = @(x) [x^2; 1+x^2];   % vanilla w/ vectorial func for f
[xm fm] = fminbnd(@(x) min(f(x)), x(1), x(3), optimset('display','iter'))
o.xtol = 1e-12; o.verb = 1; for i=1:3, y(:,i) = f(x(i)); end
[xm fm] = evp.iterparabolafit(f, x, y, o)
clear y

% need to test cases where eps is as large as the interval, etc !
f = @(x) (exp(x)-1)^2;   % skew minimum
[xm fm] = fminbnd(f, x(1), x(3), optimset('display','iter'))
o.xtol = 1; o.verb = 1; for i=1:3, y(i) = f(x(i)); end
[xm fm i] = evp.iterparabolafit(f, x, y, o)
