% generate reference set of eigenvalues and eigenfunctions, rf shape, Neu
% Barnett 1/12/14 based on example codes from SCA paper, and
% mpspack/examples/smoothdrummodesboyd.m
% reused 12/8/15
clear;
%N = 160; % ok for k<10
%N = 200; % ok for k<20
N = 450;  % ok k=40
%N = 720;   % ok for k<100
%N = 2100;   % ok for k<300
s = segment.smoothnonsym(N, 0.3, 0.2, 3); % create a closed curve
%s = segment.smoothnonsym(N, 0, 0, 3); % disc
if 0, k=20; Ns = 160:20:240; for i=1:numel(Ns), s.requadrature(Ns(i));
    min(eig(eye(Ns(i)) + 2*layerpot.D(k, s))), end, end % check N convergence
d = domain(s, 1); % create an interior domain
s.setbc(-1, 'N'); % Neumann BCs on inside
d.addlayerpot(s, 's');        % SLP representation
p = evp(d);       % create eigenvalue problem
o.modes = 1;

o.tol = 1e-6; p.solvespectrum([40 41], 'fd', o);
save rf_neu_ref40k41_N450_p_bdry
p.showmodes(struct('fmm',1)); %  for fun


%o.maxslope=1.5; o.tol = 1e-12; p.solvespectrum([10 20], 'ms', o); p.showmodes;
%save disc_neu_ref10k20_N220ms_p_bdry   % 'ms' since disc has double EVs

%p.solvespectrum([90 100], 'fd', o); save rf_neu_ref90k100_N720_p_bdry
%p.solvespectrum([300 302], 'fd', o); save rf_neu_ref300k302_N2100_p_bdry

%p.solvespectrum([10 20], 'fd', o); %p.showmodes;

