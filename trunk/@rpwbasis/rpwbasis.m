classdef rpwbasis < handle & basis

% RPWBASIS - create a real (as opposd to evanescent) plane wave basis set
%
%  b = RPWBASIS(N, k, opts) creates a real plane wave basis
%   object, with N directions, wavenumber k, and options:
%   opts.real: if true, real values (cos/sin type), otherwise complex exp.

  properties
    real     % true if sin/cos, false for complex exponentials
    dirs     % row vec of directions as C-numbers on unit circle
  end

  methods
    function b = rpwbasis(N, k, opts) % ............ constructor
      if nargin<1, N = 20; end        % Default degree of FB fct.
      if nargin<2, k = NaN; end
      if nargin<3, opts = []; end
      if ~isfield(opts,'real'), opts.real = 1; end   % default is real
      b.N = N; b.k = k; b.real = opts.real;
      b.Nf = 2*N;                     % since two funcs per direction
      b.dirs = exp(1i*pi*(1:N)/N);    % uniform in angle (need not be)
    end
    
    function [A Ax Ay] = eval(b, p, opts) % .............. evaluator at points p
      N = b.N;
      ks = b.k * b.dirs;               % kvectors as C-#s
      kdotx = real( repmat(ks, size(p.x)) .* conj(repmat(p.x, size(ks))) );
      c = cos(kdotx); s = sin(kdotx);  % NB kdotx is matrix of k dot x
      clear kdotx
      if b.real
        A = [c s];
      else
        A = c+1i*s; A = [A conj(A)];   % fwd & bkd travelling waves
      end
      if nargout==2
        kdotn = real( repmat(ks, size(p.nx)) .* conj(repmat(p.nx, size(ks))) );
        if b.real
          Ax = [-kdotn.*s kdotn.*c];   % compute An
        else
          Ax = 1i*kdotn.*A(:,1:N);
          clear kdotn
          Ax = [Ax conj(Ax)]; % fwd & bkd trav waves
        end
      elseif nargout==3
        if b.real
          kc = real( repmat(ks, size(p.nx)) );  % x-component of k vectors
          Ax = [-kc.*s kc.*c];
          kc = imag( repmat(ks, size(p.nx)) );  % y-component of k vectors
          Ay = [-kc.*s kc.*c];
        else
          kc = real( repmat(ks, size(p.nx)) );  % x-component of k vectors
          Ax = 1i*kc.*A(:,1:N); Ax = [Ax conj(Ax)];
          kc = imag( repmat(ks, size(p.nx)) );  % y-component of k vectors
          Ay = 1i*kc.*A(:,1:N);
          clear kc
          Ay = [Ay conj(Ay)];
        end
      end
    end % func
    
    function showgeom(bas, opts) % .................. crude show directions
      if nargin<2, opts = []; end
      ze = zeros(size(bas.dirs));
      plot([ze; real(bas.dirs)], [ze; imag(bas.dirs)], 'r-');
      if isfield(opts, 'label')
        n = ceil(numel(bas.dirs)/2);
        text(real(bas.dirs(n)), imag(bas.dirs(n)), opts.label);
      end
    end % func
    
  end % methods
end
