% QPUNITCELL - create a quasi-periodic unit cell domain object
%
%  uc = qpunitcell(e1, e2, k, M) creates a quasi-periodic trapezium unit cell
%   centered at the origin, and C-numbers e1 and e2 defining the
%   `horizontal' and `vertical' lattice vectors (e1 cross e2 must be positive).
%   M (which may be []) defines quadrature size (Note: number of points may be
%   M or M+1 depending on quadrature type). k is the overall wavenumber.
%
%  uc = qpunitcell(e1, e2, k, M, quad) sets quadrature type as in SEGMENT.
%
%  See also SETBLOCH, DOMAIN
classdef qpunitcell < handle & domain
  properties
    a, b                 % alpha, beta: QP phase factors
    kbloch               % Bloch wavevector (as C-#)
    e1, e2               % lattice displacement vectors (C-#s)
    r1, r2               % reciprocal lattice vectors (C-#s)
    recip                % 2x2 recip lattice matrix
    L, B, R, T           % L, B, R, T segments
  end
  
  methods
    function uc = qpunitcell(e1, e2, k, M, quad) % ................ constructor
      if nargin<4 | isempty(M), M = 20; end       % default M
      if nargin<5 | isempty(quad), quad='g'; end  % default quadrature type
      
      L = segment(M, [0 e2] - (e1+e2)/2, quad);
      B = segment(M, [e1 0] - (e1+e2)/2, quad);
      R = translate(L, e1);
      T = translate(B, e2);
      uc = uc@domain([L B R T], [-1 -1 1 1]); % create uc as domain instance
      uc.L = L; uc.B = B; uc.R = R; uc.T = T; % note LBRT normals neq domain's!
      uc.recip = 2*pi*inv([real(e1) real(e2); imag(e1) imag(e2)]);
      uc.e1 = e1; uc.e2 = e2;
      uc.r1 = uc.recip(1,1) + 1i*uc.recip(2,1);   % reciprocal lattice vecs
      uc.r2 = uc.recip(1,2) + 1i*uc.recip(2,2);
      uc.k = k; uc.setbloch;                      % default Bloch vector
    end % func
    
    function setbloch(uc, a, b)
    % SETBLOCH - choose Bloch wavevector or alpha,beta phases
    %
    %  setbloch(uc, kbloch) sets Bloch wavevector as the C-number kbloch
    %
    %  setbloch(uc, a, b) sets Bloch wavevector using the alpha, beta phase
    %   factor parameters (a,b on unit circle).
    %
    %  setbloch(uc) sets default values
      if nargin<2
        uc.kbloch = 0; uc.a = 1; uc.b = 1;     % default
      elseif nargin<3                          % interpret a as Bloch vector
        uc.kbloch = a;
        uc.a = exp(1i*real(conj(a)*uc.e1));
        uc.b = exp(1i*real(conj(a)*uc.e2));
      else                                 % find Bloch vec (in Brillouin zone)
        if abs(abs(a)-1)+abs(abs(b)-1) > 1e-15
          fprintf('setbloch warning: a or b is not on the unit circle!\n')
        end
        uc.a = a; uc.b = b;
        uc.kbloch = (uc.r1 * angle(a) + uc.r2 * angle(b))/(2*pi);
      end
    end % func
    
    function plot(uc)   % ..................... overloads domain plot
      h = domain.showsegments([uc.L uc.B uc.R uc.T], 1); % show 4 segments
    end
    
    function showbrillouin(uc)
    % SHOWBRILLOUIN - k-space plot with Brillouin zone and Bloch wavevector
      g = gcf; figure(g); hold on;
      bz = uc.recip * ([0 0 1 1 0; 0 1 1 0 0] - .5); % transformed unit square
      plot(bz(1,:), bz(2,:), 'k-'); xlabel('k_x'); ylabel('k_y');
      if ~isempty(uc.kbloch)
        hold on;
        plot([0 real(uc.kbloch)], [0 imag(uc.kbloch)], '-', 'linewidth', 5);
        plot(real(uc.kbloch), imag(uc.kbloch), 'r.', 'markersize', 30);
      end
      hold on; plot(0, 0, 'k.', 'markersize', 20);
    end % func
    
    function [Q data] = evalbasesdiscrep(uc, opts)
    % EVALBASESDISCREP - fill Q matrix mapping UC bases coeffs to discrepancy
    %
    % * to do: save data as cell array for each basis obj, reuse if opts.data!
      uc.noff = 0;                % since only one domain
      Q = [];
      for b=uc.bas                % loop over basis set objects in unit cell
        bas = b{1};               % ugly, extracts object from cell
        Qb = bas.evalunitcelldiscrep(uc);
        Q = [Q Qb];               % stack as columns in basis' order
      end
    end % func

    function addqpuclayerpots(uc, k, opts) % ...................................
    % ADDQPUCLAYERPOT - adds a sticking-out layer potential basis to unit cell
    %
    %  makes 4 instances of qpuclayerpot
      if nargin==1 | isempty(k), k = uc.k; end
%      uc.bas = {qpuclayerpot(uc, 'L', 'd', k)};
      uc.bas = {qpuclayerpot(uc, 'L', 'd', k), qpuclayerpot(uc, 'L', 's', k), qpuclayerpot(uc, 'B', 'd', k), qpuclayerpot(uc, 'B', 's', k)};
    end
    
  end % methods
end
