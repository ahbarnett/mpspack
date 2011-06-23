% Example of computing a reference set of Laplacian eigenvalues and eigenmodes
% of a smooth 2D drum using BIE and Boyd rootfinding. Dirichlet & Neumann case.
% Barnett 6/10/11.  Also see ../test/testevp.m for plotting methods

clear; a = 0.3; b = 0.2; w = 3; % shape params for generic smooth drum
N = 160;                        % # quadrature nodes good up to k=11
s = segment.radialfunc(N, {@(q) 1 + a*cos(w*(q+b*cos(q))), ...
                    @(q) -a*sin(w*(q+b*cos(q))).*w.*(1-b*sin(q)), ...
                    @(q) -a*cos(w*(q+b*cos(q))).*w^2.*(1-b*sin(q)).^2 + ...
                    a*sin(w*(q+b*cos(q))).*w.*b.*cos(q)}); % includes curvature
d = domain(s,1);                % create an interior domain
if 1 % for homog Dirichlet BCs use...
  s.setbc(-1, 'D');             % BC's applied on inside: note -1
  d.addlayerpot(s, 'd');        % DLP representation
else % for homog Neumann BCs use...
  s.setbc(-1, 'N');             % BC's applied on inside: note -1
  d.addlayerpot(s, 's');        % SLP representation
end
p = evp(d);                     % sets up eigenvalue problem object

% choose a correction scheme for the periodic quadrature...
d.bas{1}.quad = 'm';           % Kress (default, slower to fill, most accurate)
%d.bas{1}.quad = 'a'; d.bas{1}.ord = 16; % 16th-order Alpert band-diagonal

% Solve all modes in 2<k<11: (opts struct here requests mode boundary-funcs)
p.solvespectrum([2 11], [], struct('modes',1));
% the object p now contains p.kj eigenfrequencies (sqrt of eigenvalues),
% p.ndj the boundary functions (normal-derivatives for Dirichlet case),
% and error measures as follows: p.err.ej the Boyd root imaginary parts,
% p.err.minsigj the singular values of (I +- D) at the kj.

% An example higher-k calculation is 2<k<100 which requires N=720 for accuracy
% close to machine precision.
% k around 300 requires N=2500 with 16th-order Alpert.
% k around 1000 requires N=8000 "
