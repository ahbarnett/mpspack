% test segment class
% barnett 7/7/08
%
% Also see for more thorough tests: TESTDOMAIN

l = segment(10, [0 1+1i]);
figure; subplot(1,2,1); l.plot(1); axis equal; title('line seg with pm=+');
subplot(1,2,2); l.plot(-1); axis equal; title('line seg with pm=-');

a = segment(20, [1 1 pi/4 7*pi/4]);    % CCW arc
figure; subplot(1,2,1); a.plot(1); axis equal; title('CCW arc seg with pm=+');
subplot(1,2,2); a.plot(-1); axis equal; title('CCW arc seg with pm=-');

b = segment(20, [1 1 pi/4 -pi/4]);    % CW arc
figure; subplot(1,2,1); b.plot(1); axis equal; title('CW arc seg with pm=+');
subplot(1,2,2); b.plot(-1); axis equal; title('CW arc seg with pm=-');

% analytic periodic segment with periodic trapezoid quadrature...
a = 0.3; w = 3;                      % strength of variation, # wobbles.
R = @(t) 1 + a*cos(w*t); Rt = @(t) -w*a*sin(w*t);
Z = @(s) exp(2i*pi*s).*R(2*pi*s);
z = segment(100, {Z, @(s) 2*pi*(1i*Z(s) + exp(2i*pi*s).*Rt(2*pi*s))}, 'p');
figure; subplot(1,2,1); z.plot(1); axis equal; title('analytic seg with pm=+');
subplot(1,2,2); z.plot(-1); axis equal; title('analytic seg with pm=-');
