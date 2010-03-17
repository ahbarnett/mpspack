% test extrap: Richardson-type extrapolation to eval f(0). Barnett 3/15/10

clear all
fs{1} = @(x) 1./(1+x).^3;          % has vaguely-near pole singularity (x=-1)
fs{2} = @(x) exp(x);               % entire, better
fs{3} = @(x) 1+tan(pi/2 * x);      % pole at x=-1
fs{4} = @(x) exp(-(x+1).^-2);      % essential sing at x=-1 ...worse!
% (this is bad because f(0) << f(hmax) so 4 digits are lost there to roundoff)
fs{5} = @(x) sin(100*x+1);         % oscillatory (use hmax =< 0.03 for this)
fs{6} = @(x) sin((10:10:100)*x+1); % row vector test case (same as above)

hmax = 0.1;

for i=1:numel(fs)        % for each function, test 3 calling modes...
  f = fs{i}; fprintf('\n%s:\n', char(f))
  disp('default...');
  [f0 err N] = utils.extrap(f, hmax);
  fprintf('actual err = %.3g, pred err = %.3g, N = %d.\n', ...
          max(abs(f0-f(0))), max(err), N)
  disp('fix reltol=1e-6...');
  [f0 err N] = utils.extrap(f, hmax, struct('reltol', 1e-6));
  fprintf('actual err = %.3g, pred err = %.3g, N = %d.\n', ...
          max(abs(f0-f(0))), max(err), N)
  disp('fix N=10...');
  [f0 err N] = utils.extrap(f, hmax, struct('N', 10));
  fprintf('actual err = %.3g, pred err = %.3g, N = %d.\n', ...
          max(abs(f0-f(0))), max(err), N)
end
