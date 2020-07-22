% test smoothfourierz
% Barnett 7/22/20
clear

% band-limited, should be exact
n = 10;
s0 = segment.smoothstar(n,0.3,3);
z = s0.x;
s = segment.smoothfourierz(100,z);
figure; plot(s0.x,'.','markersize',20); hold on; s.plot;
s1 = segment.smoothstar(100,0.3,3);


% lens, not C^1. all terms. Induces Gibbs ringing
n = 40;
z = exp(2i*pi*(-n/4:n/4-1)/(2*n)) - sqrt(.5);
z = [z -z];
s = segment.smoothfourierz(200,z);
figure; s.plot; hold on; plot(z,'.','markersize',20); title('lens, Gibbs')

% smooth roll-off of Fourier coeffs, localizes the ringing
s = segment.smoothfourierz(200,z,n/2,1e-6);
figure; s.plot; hold on; plot(z,'.','markersize',20); title('lens, rolloff');
% haven't tested how accurate the spatial localization is, ie how well
% fits distant points in locally high-order-smooth region.


