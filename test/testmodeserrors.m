% test the evp.modeserrors function. Barnett 12/21/11
clear all classes
kwin = [30 31]; N = 300;
s = segment.smoothstar(N, 0.3, 5); % pentafoil
d = domain(s, 1);                % create an interior domain
s.setbc(-1, 'D');               % Dirichlet BC's applied on inside: note -1
p = evp(d);                  % sets up eigenvalue problem object
q = utils.copy(p);

o.modes = 1; o.khat = 'r'; o.fhat = 's';  % compare two NtD-scaling method runs
o.eps = 0.1; p.solvespectrum(kwin, 'ntd', o);
o.eps = 0.2; q.solvespectrum(kwin, 'ntd', o);

p.kj-q.kj                              % compare wavenumbers
o = []; o.wei = 1; p.modeserrors(q,o); % compare modes (some zero, some 1e-3)
