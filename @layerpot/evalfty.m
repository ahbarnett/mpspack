function [A Ax] = evalfty(b, f, o) % ...........for 1.5D scatt, qpstrip
% EVALFTY - evaluate layer potential on y-Fourier Sommerfeld contour nodes
%
% A = EVALFTY(b, f) gives matrix A mapping layer potential density
%   values to values at the complex k_y points for the vertical y-slice
%   defined by ftylayerpot basis f.
%
% [A Ax] = EVALFTY(b, f) also returns x-derivatives (normal derivs on
%   the y-slice).
%
% [A Ax] = EVALFTY(b, f, opts) allows certain options, such as:
%   opts.reflect = 'd' or 'i' adds in effect of an image reflected through y=h
%   opts.reflecth = h value to reflect about (required if opts.reflect used).
%   opts.reflecteta = eta in u_n+i.eta.u=0, for opts.reflect='i' (default om).
%
% Notes: 
% 1) b is only used to extract k_y and origin information, not to
%    evaluate as a basis set.
% 2) No close-evaluation is used - in fact there is no quadrature error here.
% 3) Signs in the formulae were guessed by trial and error.
%
% Part of the 1.5D QP scatt tools.
%
% See also: TESTLAYERPOTEVALFTY

% (C) Alex Barnett 2010 - 2011
if nargin<3, o = []; end
if ~isfield(o,'reflect'), o.reflect = ''; end % default no reflection (false)
if o.reflect, if ~isfield(o,'reflecth'), error('needs opts.reflecth!'); else
  h = o.reflecth; end, end
om = b.k(o);                                   % wavenumber (from layerpot)
d = b.seg.x - f.orig;                          % displacements from slice origin
x = real(d); y = imag(d);                      % Cartesian displ, col vecs
nx = real(b.seg.nx); ny = imag(b.seg.nx);      % x,y cpts of source normal
N = numel(d);                                  % # src pts
k = f.kj;                                      % row vec
M = numel(k);                                  % # target complex k_y wave#s
som = sqrt(om^2 - k.^2);                       % Sommerfeld, row vec
exf = exp(1i*som.'*abs(x).');                  % decay matrix (outer prod)
pha = exp(-1i*k.'*y.');                        % y-translation mat (outer prod)

B = pha .* exf .* repmat(b.seg.w/4/pi, [M 1]); % common to all cases (inc 1/4pi)

% the following could be sped up my testing if b.a(1)==0 or b.a(2)==0 ...
A = B .* (repmat(b.a(1)*1i./som.', [1 N]) + ...        % SLP
          repmat(-b.a(2)*nx.'.*sign(x).', [M 1]) + ... % DLP x-deriv part
          (b.a(2)*k./som).'*ny.');                     % y-deriv via outer prod
if nargout>1
  Ax = B .* (repmat(b.a(1)*sign(x).', [M 1]) + ...     % SLP
             (b.a(2)*1i*som.')*nx.' - ...              % DLP x-deriv, outer prod
             (b.a(2)*1i*k.')*(ny.*sign(x)).');         % y-deriv via outer prod
end
if o.reflect                           % Half-space reflection cases:
  pha = exp(-1i*k.'*(2*h-y.'));        % (reflect off y=h, same above or below)
  if isfield(o,'reflecteta'), eta = o.reflecteta; else eta = om; end % pick eta
  if o.reflect=='d'                 % add in for Dirichlet half-plane Green's
    B = -pha .* exf .* repmat(b.seg.w/4/pi, [M 1]);  % note overall minus for d
  elseif o.reflect=='i'                % add in impedance half-plane Green's
    B = pha .* exf .* repmat(b.seg.w/4/pi, [M 1]) ...
        .* repmat(((k-eta)./(k+eta)).', [1 N]);   % the impedance case
  end
  A = A + B .* (repmat(b.a(1)*1i./som.', [1 N]) + ...        % SLP
              repmat(-b.a(2)*nx.'.*sign(x).', [M 1]) - ...   % DLP x-deriv part
              (b.a(2)*k./som).'*ny.');                       % y-deriv
  if nargout>1
    Ax = Ax + B .* (repmat(b.a(1)*sign(x).', [M 1]) + ...     % SLP
             (b.a(2)*1i*som.')*nx.' + ...              % DLP x-deriv, outer prod
             (b.a(2)*1i*k.')*(ny.*sign(x)).');         % y-deriv via outer prod
  end
end
