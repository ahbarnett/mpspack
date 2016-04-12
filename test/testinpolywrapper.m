% test code for inpolygon wrapper routine, validating against matlab's native.
% This is OBSOLETE since switched to MATLAB's native inpoly in around 2013.
% It is retained for amusement only.

% Copyright (C) 2008 - 2012, Timo Betcke, Alex Barnett

wrapperusesmatlab = 0;             % flag, set to 0 if you edit inpolywrapper
                                   % to use much faster MEX file.
% pathological examples
i = utils.inpolywrapper([], [0; 1; 1+1i; 1i])
if isempty(i), disp('ok'); else warning('wrong'); end
if ~wrapperusesmatlab
  i = utils.inpolywrapper([0.5+0.5i; -.5+.5i], [])  % fails for matlab inpolygon
  if i==[0;0], disp('ok'); else warning('wrong'); end
  i = utils.inpolywrapper([0], [])           % test point is pure real
  if i==0, disp('ok'); else warning('wrong'); end
end
i = utils.inpolywrapper([0], [0])          % test point and polygon pure real
if i==[0], disp('ok'); else warning('wrong'); end
if ~wrapperusesmatlab
  i = utils.inpolywrapper([], [])
  if isempty(i), disp('ok'); else warning('wrong'); end
end

% small example
i = utils.inpolywrapper([0.5+0.5i; -.5+.5i], [0; 1; 1+1i; 1i])
if i==[1;0], disp('ok'); else warning('wrong'); end

% medium example w/ multiple calls
nv = 100;
N = 1e3;
n = 1e2; % # of calls
v = rand(nv,1).*exp(2i*pi*(1:nv)'/nv);    % random nv-gon 
p = rand(N,1) + 1i*rand(N,1);  % all inside bounding box
disp('medium tests (100% lie in bounding box): matlab...')
tic; for j=1:n, ii = inpolygon(real(p),imag(p), real(v), imag(v)); end, toc
fprintf('vertices-points/sec = %.3g\n',nv*N*n/toc)
disp('inpolywrapper...');
tic; for j=1:n, i = utils.inpolywrapper(p, v); end, toc
fprintf('vertices-points/sec = %.3g\n',nv*N*n/toc)
nwrong = numel(find(i~=ii));
if nwrong==0, disp('ok'); else warning(sprintf('nwrong=%d\n',nwrong)); end
if wrapperusesmatlab, disp('above timings will be the same, redundant'); end

% big example (timings on i7-3720QM, w/ bb test in inpolywrapper):
% matlab 7-9s, Engwirda 0.33s, Franklin 0.2s, Luong 0.08s
nv = 1000;
N = 1e6;
v = rand(nv,1).*exp(2i*pi*(1:nv)'/nv);    % random nv-gon 
p = sqrt(10)*rand(N,1) + 1i*sqrt(10)*rand(N,1); % note most outside bounding box
disp('big test (10% lie in bounding box): matlab...')
tic; ii = inpolygon(real(p),imag(p), real(v), imag(v)); toc
fprintf('vertices-points/sec = %.3g\n',nv*N/toc)
disp('inpolywrapper...');
tic; i = utils.inpolywrapper(p, v); toc
fprintf('vertices-points/sec = %.3g\n',nv*N/toc)
nwrong = numel(find(i~=ii));
if nwrong==0, disp('ok'); else warning(sprintf('nwrong=%d\n',nwrong)); end
if wrapperusesmatlab, disp('above timings will be the same, redundant'); end

% another big eg, all in bb
nv = 1000;
N = 1e5;
v = rand(nv,1).*exp(2i*pi*(1:nv)'/nv);    % random nv-gon 
p = 2*rand(N,1) + 2i*rand(N,1) - 1 - 1i;  % note all in bb
disp('big test (100% lie in bounding box): matlab...')
tic; ii = inpolygon(real(p),imag(p), real(v), imag(v)); toc
fprintf('vertices-points/sec = %.3g\n',nv*N/toc)
disp('inpolywrapper...');
tic; i = utils.inpolywrapper(p, v); toc
fprintf('vertices-points/sec = %.3g\n',nv*N/toc)
nwrong = numel(find(i~=ii));
if nwrong==0, disp('ok'); else warning(sprintf('nwrong=%d\n',nwrong)); end
if wrapperusesmatlab, disp('above timings will be the same, redundant'); end
