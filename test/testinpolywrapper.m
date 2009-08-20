% test code for inpolygon wrapper routine, validating against matlab's native.
% Copyright (C) 2008, 2009, Timo Betcke, Alex Barnett

wrapperusesmatlab = 1;             % flag, set to 0 if you edit inpolywrapper
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

% big example
nv = 200;
N = 1e5;
v = exp(2i*pi*(1:nv)'/nv);    % nv-gon on unit circle
p = rand(N,1) + 1i*rand(N,1);
disp('matlab...')
tic; ii = inpolygon(real(p),imag(p), real(v), imag(v)); toc
disp('inpolywrapper...');
tic; i = utils.inpolywrapper(p, v); toc
nwrong = numel(find(i~=ii));
if nwrong==0, disp('ok'); else warning(fprintf('nwrong=%d\n',nwrong)); end
if wrapperusesmatlab, disp('above timings will be the same, redundant'); end
