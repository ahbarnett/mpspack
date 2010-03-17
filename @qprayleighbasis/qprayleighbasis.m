classdef qprayleighbasis < handle & basis

% QPRAYLEIGHBASIS - create a radiative Rayleigh Fourier series in half-strip
%
% b = QPRAYLEIGHBASIS(seg, pm, N) creates a Rayleigh series basis with radiation
%   condition facing opposite to (if pm=+1) or along (if pm=-1) segment normal,
%   based upon the horizontal segment seg. N is the order of the basis
%   (there are 2N+1 functions).

% (C) 2010 Alex Barnett

  properties
    yo                               % y-value where unit-normalized
    Lo, Ro                           % x-values of left and right sides
    seg, pm                          % segment and sense (as in qpstrip)
    up                               % sense of radiation condition
    d                                % width = Ro-Lo
    a                                % x-Bloch phase alpha
  end
  
  methods
    function r = qprayleighbasis(s, pm, N, opts)   % constructor
      if nargin<4, opts=[]; end
      if nargin<3 || isempty(N), N=20; end                  % default order
      if nargin<2, up = +1; end
      if nargin<1, s = segment([], [0 1]); end              % default seg
      if abs(imag(s.eloc(1)-s.eloc(2)))>1e-14
        error('seg must be horizontal!'); end
      r.yo = imag(s.eloc(1));
      r.Lo = min(real(s.eloc)); r.Ro = max(real(s.eloc)); r.d = r.Ro-r.Lo;
      if real(s.eloc(1))<real(s.eloc(2))        % L-to-R seg
        r.up = pm;
      else, r.up = -pm; end
      r.N = N; r.seg = s; r.pm = pm; r.nmultiplier = 1;     % default nmult
    end
    
    function h = showgeom(r, o)
      h = plot(r.seg, r.pm, struct('arrow',0)); set(h,'color', [1 0 0]);
    end
    
    function Nf = Nf(r) % ...................... returns # dofs
      Nf = 2*r.N + 1;
    end
        
    function [A, A1, A2] = eval(r,p,o) % ........................ evaluator
        % EVAL - evaluates Rayleigh Fourier basis at given set of points
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
        %    derivatives in the x- and y-directions. That is, Ax has ij'th
        %    entry d/dx Phi_j(z_i) while Ay has ij'th entry d/dy Phi_j(z_i)
        %
        % Also see: POINTSET, QPSTRIP
      N = r.N; ns = -N:N;            % Fourier series indices (row vec)
      M = length(p.x);               % Number of eval points
      x = real(p.x) - r.Lo;          % x-displacements from left end origin
      y = imag(p.x) - r.yo;          % y-displacements from origin
      om = r.k;                      % inherits wavenumber from domain
      kx = (real(log(r.a)/1i) + 2*pi*ns) / r.d; % x-wavenumbers (shift by Bloch)
      ky = r.up * sqrt(om^2 - kx.^2);     % y-wavenumbers (w/o i, also row vec)
      A = exp(1i*(x*kx + y*ky));     % outer prods to eval
      A = A .* repmat(1./max(abs(ky),1), [M 1]);     % rescale using deriv
      if nargout==2
        A1 = A .* ((1i*real(p.nx))*kx + ((1i*imag(p.nx))*ky)); % n-deriv
      elseif nargout==3
        A1 = A .* repmat(1i*kx, [M 1]);                        % x-,y-derivs
        A2 = A .* repmat(1i*ky, [M 1]);
      end
    end
    
  end % ... end methods
end   
