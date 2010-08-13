function [x derr y u ier] = intervalrootsboyd(f, int, o)
% INTERVALROOTSBOYD - find roots of real analytic function to spectral accuracy
%
% x = intervalrootsboyd(f, [a b]) finds roots of the function f in the real
%   interval (a b], using sampling at Chebychev points which are doubled in
%   resolution until Fourier coefficients have sufficiently decayed, then via
%   J. Boyd's 2002 method of degree doubling. Method is O(N^3) in the number of
%   samples used, due to the use of Matlab roots command. Idea from CHEBFUN
%   package of Trefethen, Driscoll, Platte, et al.
%   Notes:
%   f need not accept vector of input arguments - it is only called as scalar.
%   Zeros precisely at endpoints are not guaranteed to be included or excluded.
%
% x = intervalrootsboyd(f, [a b], opts) allows certain options such as:
%   opts.Ftol  : Fourier coeff decay relative tolerance (default 1e-12)
%   opts.Nmax : largest number of function evaluations (default 256)
%   opts.disp : 0,1,... verbosity
% Other opts are passed to trigpolyzeros
%

% [x derr y u ier] = intervalrootsboyd(f, [a b], ...) also returns estimated
%   errors in the form of distances of roots from the real axis, the final set
%   of ordinates y, and function values u used, and error flag (0 if success,
%   otherwise gives failure mode).
%
% Also see: TRIGPOLYZEROS

% (C) 2010, Alex Barnett

if nargin<3, o = []; end
if ~isfield(o, 'Ftol'), o.Ftol = 1e-12; end
if numel(o.Ftol)~=1, error('opts.Ftol must be a single number!'); end
o.real = 1;            % for trigpolyzeros, complex zeros but only [0,pi)
if ~isfield(o, 'Nmax'), o.Nmax = 256; end  % means roots does a 1024-sized eig
if ~isfield(o, 'disp'), o.disp = 0; end
a = int(1); b = int(2); rad = (b-a)/2; cen = (b+a)/2;

N = 4;                 % half minimum number of points on half-circle
t = pi*(0:N)/N;        % angles theta
u = nan(size(t));      % the data
y = cen + rad*cos(t);  % ordinates
for i=1:numel(y), u(i) = f(y(i)); if o.disp, fprintf('f(%.16g)=%.16g\n',y(i),u(i)), end, end  % initialize func evals
F = [1 1];                 % dummy value greater than tol
Fmetr = @(F) max(abs(F(1:2)))/max(abs(F));  % metric for decay to small F coeffs
x = []; derr = []; 

while N<o.Nmax && Fmetr(F)>o.Ftol  % large coeffs small?
  N = 2*N; t = pi*(0:N)/N; y = cen + rad*cos(t); u(1:2:N+1) = u; % reuse f-evals
  for i=1:N/2, u(2*i) = f(y(2*i)); if o.disp, fprintf('f(%.16g)=%.16g\n',y(2*i),u(2*i)), end, end     % fill in missing evals
  ufold = [u(end:-1:2) u(1:end-1)];        % samples on -pi+2*pi*(0:N-1)/N
  %figure; plot(ufold); drawnow  % for debugging
  F = fftshift(fft(fftshift(ufold)));         % compute shifted Fourier coeffs
  if Fmetr(F)<=o.Ftol
    [tz derr] = utils.trigpolyzeros(F, o);      % get theta roots in [-pi,pi)
    x = cen + rad*cos(tz);
  end
end
%figure; plot(abs(F)); title('|F_j|');   % debug

ier = 0;                                 % output info
if Fmetr(F)>o.Ftol
  warning('stopped even though Fourier coeffs not sufficiently decayed!');
  ier = 1;
end