% Examples of computing Laplacian Dirichlet eigenvalues and eigenmodes
% of a smooth 2D drum using the weighted-Neumann-to-Dirichlet scaling method
% published by Barnett-Hassell CPAM 2014 (see this paper for algorithm).
% Barnett 6/10/11, cleaned up 4/13/16

clear all classes; N = 300;               % # quadrature nodes good up to k=30
s = segment.smoothnonsym(N, 0.3, 0.2, 3); % closed smooth non-symm segment
d = domain(s, 1);                % create an interior domain
s.setbc(-1, 'D');               % Dirichlet BC's applied on inside: note -1
p = evp(d);                     % sets up eigenvalue problem object
p.set

% clean up this file by moving to expt_*.m 

if 0   % check operator spectrum and timing at a single k:
  tic; l = p.NtDspectrum(30); toc, min(l) % test
  tic; ll = p.NtDspectrum(30, struct('quad','a','ord',16)); toc, min(ll)
  tic; ll = p.NtDspectrum(30, struct('quad','m')); toc, min(ll)
end
  
%o.khat='l'; o.fhat='f'; p.solvespectrum([30 31], 'ntd', o); % linear, worse
%o.khat='o'; o.fhat='f'; p.solvespectrum([30 31], 'ntd', o); % 2nd order, best

% solve eigenvalues only
%p.solvespectrum([30 31], 'ntd', struct('eps',0.1,'khat','o','fhat','s'));

% solve with modes & show them
p.solvespectrum([30 31], 'ntd', struct('eps',0.1,'khat','o','fhat','s','modes',1));
p.kj  % eigenwavenumbers found
p.showmodes;


if 0 % try higher k range now, medium-scale calc......
  s.requadrature(720);
  tic; p.solvespectrum([90 100], 'ntd', struct('eps',0.1,'khat','l')); toc;
  % 1271 sec for eigfreqs only (21 mins)
  % 12.25 sec per kstar eval, ie 2.5 sec per eigfreq found.
  % 1804 sec with eigfuncs (30 mins, ie 1.5 times longer).
  tic; p.solvespectrum([90 100], 'ntd', struct('eps',0.1,'khat','l','modes',1,'fhat','f')); toc;
  %save solvespec_90k100_khatl_fhatf.mat p
end

if 0   % timings of aspects and other methods at this medium k.....
  d.addlayerpot(s, 'd');        % basis for Dir
  tic;  p.solvespectrum([90 90.1], 'fd'); toc;
  %ratio 5 evals per mode found
  % 428 sec for 90 to 91, so 4300 sec = 72 mins, just eigvals
  % 1.65 sec per fred det eval, 8.7 sec per eigval found.
  tic;  p.solvespectrum([90 100], 'fd', struct('modes',1)); toc;
  % If want eigvecs and use dense svd, extra 15 sec per eigvec found.
  % If want eigvecs and use iter=1 w/ lu, only extra 1.8 sec per eigvec found.
  tic; ll = p.NtDspectrum(90, struct('quad','m')); toc % 12 sec
  % dominated by 8.5 sec for complex eig (vs 0.4 for complex det or lu).
  % 14.5 s if want eigvecs too. complex inv is only 1.2 sec, complex svd 3.1 sec
  % or more like 15 if want U and V.
  profile clear; profile on; ll = p.NtDspectrum(90, struct('quad','m')); profile off; profile viewer
end
