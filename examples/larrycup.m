function s = larrycup(N)
% parametrized resonant cup analytic curve. Barnett 9/14/14
% for mpspack format, via Fourier series,  3/6/19.

a = 0.2;  % thickness: cannot be >1/3 otherwise not smooth to emach
b = pi/6;  % controls approx opening angle in radians (keep small for resonant)

N = ceil(N/2);
s = ((1:N)-0.5)/N * pi;  % note half-offset, needed for easy reflection abt z
r = 1 - a*erf((s-pi/2)/a);  % radius: starts at 1+a, ends at 1-a
c = a; %*(1-b/pi);  % is theta rounding scale
sabs = @(x) exp(-(x/c).^2)*c/sqrt(pi)+x.*erf(x/c); % c-smoothed absval
th = b-a + 2*(1-(b-a)/pi)*sabs(s-pi/2);
%th = b + sqrt(log(exp((2*(1-b/pi)*(s-pi/2)).^2) + exp(c^2))); % theta

% coords in (rho,z) plane:
rho = r.*sin(th); z = r.*cos(th);  % theta down from z axis as in 3D cyl coords

z = z*1.2;  % vert stretch! makes ellipse cavity

s = segment();
Z = [rho -rho(end:-1:1)] + 1i*[z z(end:-1:1)]; % complex coords of full curve
N = numel(Z);
% (appropriate for half-integer offset
%figure; semilogy(abs(fft(Z))); title('Fourier coeff decay, to close to emach?')
%Z = Z(end:-1:1);
zhat = fft(Z(:))/N;

s = segment(N,{@(t) fourierZ(zhat,t), @(t) fourierZp(zhat,t), @(t) fourierZpp(zhat,t)},'p');

%figure; plot(Z,'k.'); hold on; l=0.1; plot([Z;Z+l*Zn],'b-'); % show it
%axis equal xy tight;



% analytic formulae for a Fourier segment --------------

function z = fourierZ(zhat,t)     % must work on vector of t's
t = 2*pi*t;
N =numel(zhat);  % even
z = 0*t;
for k=0:N/2
  z = z + zhat(k+1)*exp(1i*k*t);
end
for k=-N/2:-1
  z = z + zhat(k+1+N)*exp(1i*k*t);
end

function zp = fourierZp(zhat,t);  % deriv func Z'
N = numel(zhat);
zp = 2*pi*fourierZ(zhat.*[0 1i*(1:N/2-1) 0 1i*(-N/2+1:-1)].', t);

function zpp = fourierZpp(zhat,t);  % deriv func Z''
N = numel(zhat);
zpp = 2*pi*fourierZp(zhat.*[0 1i*(1:N/2-1) 0 1i*(-N/2+1:-1)].', t);

% ---------------------

function g = perispecdiff(f)
% PERISPECDIFF - use FFT to take periodic spectral differentiation of vector
%
% g = PERISPECDIFF(f) returns g the derivative of the spectral interpolant
%  of f, which is assumed to be the values of a smooth 2pi-periodic function
%  at the N gridpoints 2.pi.j/N, for j=1,..,N (or any translation of such
%  points).
%
% Barnett 2/18/14
N = numel(f);
if mod(N,2)==0   % even
  g = ifft(fft(f(:)).*[0 1i*(1:N/2-1) 0 1i*(-N/2+1:-1)].');
else
  g = ifft(fft(f(:)).*[0 1i*(1:(N-1)/2) 1i*((1-N)/2:-1)].');
end
g = reshape(g,size(f));
