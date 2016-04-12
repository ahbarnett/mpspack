% Use MPSpack & Fredholm-det method for low Dirichlet eigenmodes.
% Barnett 10/5/12 for Peter Simon.

clear; verb = 0;  % set verb = 1 for output plots
bc = 'D'; % boundary conditions: 'D' or 'N' - D gets 8 digits, N 4-5 digits
o.kressq=5; %<9, corner-packing parameter for Kress reparametrization

%seg=segment.smoothstar(50,0.3,3); % warm-up shapes
%seg=segment.polyseglist(50, [1, exp(3i*pi/8), exp(5i*pi/4)], 'pc', o); % triang
%seg = segment.polyseglist(40, [am/2, am/2+1i*bm, -am/2+1i*bm, -am/2+.5], 'pc', o); % trapezoid

n = 48; % convergence param (discretization pts per segment; must be 4*integer)
am = 2; bm = 1; s = 0.5; R = 0.1; g = 0.6; % Simon's ridge waveguide params
h = 1i*(bm-g-R); % y-shift of centers
ss = segment.polyseglist(n, [s/2+h, s/2, am/2, am/2+1i*bm, -am/2+1i*bm, -am/2, -s/2, -s/2+h], 'pc', o);
seg = [ss(1:end-1) segment(n/2,[-s/2+R+h,R,pi,pi/2],'pc',o) segment(n,[-s/2+R+h+1i*R, s/2-R+h+1i*R],'pc',o) segment(n/2,[s/2-R+h,R,pi/2,0],'pc',o)];

d = domain(seg, 1);
if verb, figure; d.plot; print -depsc2 ../doc/figs/ridgeguidegeom.eps; end
seg.setbc(-1, bc);               % BC's applied on inside: note -1
p = evp(d);                     % sets up eigenvalue problem object
if bc=='D', d.addlayerpot(seg, 'd'); else, d.addlayerpot(seg, 's'); end
o.modes = 1;
o.iter = 1;  % can set to 1 if no degeneracies
o.tol = 1e-4; % make sure to capture even bad ones
kint = [1 8]; tic; p.solvespectrum(kint, 'fd', o); toc
tic; [u gx gy] = p.showmodes; toc

if verb, if bc=='D', print -depsc2 ../doc/figs/ridgeguideDlow.eps
         else, print -depsc2 ../doc/figs/ridgeguideNlow.eps
  end, end
