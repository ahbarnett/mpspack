classdef epwbasis < handle & basis

% EPWBASIS - create an evanescent plane wave basis set
%
% dirs DIFFERS SLIGHTLY FROM RPWBASIS ! HACK FOR NOW for POLESFTYLP
%
%
% TODO:
%  b = EPWBASIS(N) creates a real plane wave basis object, with N directions
%   theta_j equally spaced in (0,pi], i.e. theta_j = pi.j/N. As with all basis
%   types, the wavenumber k is determined by that of the affected domain.
%   Basis functions are:
%      {cos(k n_j.z)} j=1,..,N and {sin(k n_j.z)} j=1,..,N       (real case)
%   or {exp(i k n_j.z)} j=1,..,N and {exp(-i k n_j.z)} j=1,..,N  (complex case)
%     where the direction unit vectors are n_j = (cos theta_j, sin theta_j)
%
%   If k=0 a warning is produced since this basis does not have a useful
%   limiting form as k tends to zero (regfbbasis should be used instead).
%
% See also: DOMAIN/ADDEPWBASIS

% Copyright (C) 2010 Alex Barnett, Timo Betcke

  properties
    real     % true if sin/cos, false for complex exponentials (fwd/bck)
    half     % if non-negative, 0 (1) selects only the first (second) type
    dirs     % row vec of directions as complex numbers on unit circle
  end

  methods
    function b = epwbasis(N, opts) % ............ constructor
      if nargin<1, N = 20; end        % Default degree
      if nargin<2, opts = []; end
      if ~isfield(opts,'real'), opts.real = 1; end   % default is real
      if ~isfield(opts,'half'), opts.half = -1; end  % default is both halves
      if ~isfield(opts,'nmultiplier'), opts.nmultiplier=1; end
      b.nmultiplier=opts.nmultiplier;
      b.real = opts.real; b.half = opts.half;
      b.updateN(N);           % sets up overall N, and PW angles
    end
 
    function updateN(b,N)
        % UPDATEN - Update degree N of real plane wave basis function set
        b.N=ceil(b.nmultiplier*N);
        b.dirs = exp(1i*pi*(1:b.N)/b.N);    % uniform in angle, for now
    end
        
    function Nf = Nf(b)
      Nf = 2*b.N;                     % since two funcs per direction
      if b.half>=0, Nf = Nf/2; end    % unless use only half the basis set
    end
    
    function [A Ax Ay] = eval(b, p, opts) % .............. evaluator at points p
    % EVAL - evaluates real plane wave basis set at given set of points
    %
    % A = EVAL(p) where p is a pointset object containing M points, returns
    %   a M-by-Nf matrix whose ij'th entry is Phi_j(z_i), where Phi_j is
    %   the jth basis function, and z_i the ith point. Nf is the number
    %   of degrees of freedom in the basis set object.
    %        
    % [A An] = EVAL(p) also returns matrix An whose ij'th entry is
    %   d/dn_i Phi_j(z_i), the derivative in the ith normal direction in
    %   the pointset.
    %
    % [A Ax Ay] = EVAL(p) returns A as above, and matrices of basis function
    %   derivatives in the x- and y-directions. That is, Ax has ij'th
    %   entry d/dx Phi_j(z_i) while Ay has ij'th entry d/dy Phi_j(z_i)
    %
    % Also see: POINTSET    
      N = b.N;
      
      % get x y wavevector compts:
      w = b.doms.k;
      kx = 0*b.dirs; ky = kx;
      for i=1:numel(b.dirs), d = abs(b.dirs(i));
        if abs(abs(d-1))>1e-15, %fprintf('evan: d-1=%g\n',d-1)
          ky(i) = 1i*w*(1-d^2)/2/d; kx(i) = sqrt(w^2-ky(i)^2); % choose updown
          th = angle(b.dirs(i));
          K = [cos(th) -sin(th); sin(th) cos(th)] * [kx; ky];
          kx = K(1,:); ky = K(2,:);
        else kx(i) = w*real(b.dirs(i)); ky(i) = w*imag(b.dirs(i)); end
      end

      kdotx = kron(kx,real(p.x)) + kron(ky,imag(p.x));
      c = cos(kdotx); s = sin(kdotx);  % NB kdotx is matrix of k dot x
      clear kdotx
      h = b.half; if b.real && h>=0, error 'opts.half>=0 not yet implemented for real case!'; end
      if b.real
        A = [c s];
      else
        A = c+1i*s; if h<0, A = [A conj(A)];   % fwd & bkwd travelling waves
          elseif h>0, A = conj(A); end         % only bkwd
      end
      if nargout==2
        kdotn = kron(kx,real(p.nx)) + kron(ky,imag(p.nx));
        if b.real
          Ax = [-kdotn.*s kdotn.*c];   % compute An
        else
          Ax = 1i*kdotn.*A(:,1:N);
          clear kdotn
          if h<0, Ax = [Ax conj(Ax)];          % fwd & bkd trav waves
            elseif h>0, Ax = conj(Ax); end     % only bkwd
        end
      elseif nargout==3
        if b.real
          kc = repmat(kx, size(p.nx));  % x-component of k vectors
          Ax = [-kc.*s kc.*c];
          kc = repmat(ky, size(p.nx));  % y-component of k vectors
          Ay = [-kc.*s kc.*c];
        else
          kc = repmat(kx, size(p.nx));  % x-component of k vectors
          Ax = 1i*kc.*A(:,1:N);
          if h<0, Ax = [Ax conj(Ax)];          % fwd & bkd trav waves
            elseif h>0, Ax = conj(Ax); end     % only bkwd
          kc = repmat(ky, size(p.nx));  % y-component of k vectors
          Ay = 1i*kc.*A(:,1:N);
          clear kc
          if h<0, Ay = [Ay conj(Ay)];          % fwd & bkd trav waves
            elseif h>0, Ay = conj(Ay); end     % only bkwd
        end
      end
    end % func
    
    function showgeom(b, opts) % .................. crude show directions
    % SHOWGEOM - plot real plane wave basis set geometry information
      if nargin<2, opts = []; end
      ze = zeros(size(b.dirs));
      plot([ze; real(b.dirs)], [ze; imag(b.dirs)], 'r-');
      if isfield(opts, 'label')
        n = ceil(numel(b.dirs)/2);
        text(real(b.dirs(n)), imag(b.dirs(n)), opts.label);
      end
    end % func
    
  end % methods
end
