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
    N                    % total # degrees of freedom in QP basis sets
    basnoff              % dof offsets for QP basis sets (numeric list)
    Qpolydata            % matrix polynomial storage for each bas (*cell array)
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
      uc.recip = 2*pi*inv([real(e1) real(e2); imag(e1) imag(e2)]'); % cols
      uc.e1 = e1; uc.e2 = e2;
      uc.r1 = uc.recip(1,1) + 1i*uc.recip(2,1);   % reciprocal lattice vecs
      uc.r2 = uc.recip(1,2) + 1i*uc.recip(2,2);
      uc.k = k; uc.setbloch;                      % default Bloch vector
    end % func
    
    function [N noff] = setupbasisdofs(uc)
    % SETUPBASISDOFS - set up indices of basis degrees of freedom in unit cell
    %
    %  Similar to problem.setupbasisdofs except internal to unit cell bases
      if isempty(uc.bas), fprintf('warning: no basis sets in unit cell!\n'); end
      noff = zeros(1, numel(uc.bas));
      N = 0;                         % N will be total # dofs (cols of B)
      for i=1:numel(uc.bas)
        uc.bas{i}.updateNf;
        noff(i) = N; N = N + uc.bas{i}.Nf;    % set up dof offsets of bases
      end
      uc.N = N; uc.basnoff = noff;   % store stuff as unit cell properties
    end
    
    function s = discrepsqrtwei(uc)
    % DISCREPSQRTWEI - return sqrt of unit cell discrepancy quadrature weights
    %    
    % order matches eval.evalunitcelldiscrep...
      s = sqrt([uc.L.w, uc.L.w, uc.B.w, uc.B.w]);     % row vec
    end
    
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
    
    function h = plot(uc)   % ..................... overloads domain plot
      h = domain.showsegments([uc.L uc.B uc.R uc.T], 1); % show 4 segments
    end
    
    function showbrillouin(uc)
    % SHOWBRILLOUIN - k-space plot with UC Brillouin zone and Bloch wavevector
      g = gcf; figure(g); hold on;
      bz = uc.recip * ([0 0 1 1 0; 0 1 1 0 0] - .5); % transformed unit square
      plot(bz(1,:), bz(2,:), 'k-'); xlabel('k_x'); ylabel('k_y');
      if ~isempty(uc.kbloch)
        hold on;
        plot([0 real(uc.kbloch)], [0 imag(uc.kbloch)], '-', 'linewidth', 5);
        plot(real(uc.kbloch), imag(uc.kbloch), 'r.', 'markersize', 30);
      end
      hold on; plot(0, 0, 'k.', 'markersize', 20);
      % label the standard zone points...
      text(0,0,'\Gamma', 'fontsize', 12, 'Fontweight', 'demi');
      text(real(uc.r1)/2, imag(uc.r1)/2, 'X', 'fontsize', 12, 'Fontweight', 'demi');
      text(real(uc.r1+uc.r2)/2, imag(uc.r1+uc.r2)/2, 'M', 'fontsize', 12, 'Fontweight', 'demi');
    end % func
    
    function Q = evalbasesdiscrep(uc, opts) % .................... Q
    % EVALBASESDISCREP - fill Q matrix mapping UC bases coeffs to discrepancy
    %
    %  This stacks matrices from basis.evalunitcelldiscrep, and handles poly
    %   matrix data storage within the cell array uc.Qpolydata.
    %
    %  Q = EVALBASESDISCREP(uc) returns matrix of basis function evaluations
    %   given basis set dofs as stored in unit cell.
    %
    %  Q = EVALBASESDISCREP(uc, opts) passes in options including:
    %   opts.poly: if true, returns 3d-array Q(:,:,1:5) of polynomial coeff
    %    matrices for alpha^{-2}, alpha^{-1}, ... , alpha^2.
    %
    % To do: change block stacking to preallocation for speed, as in func below
    Q = [];
      opts.nei = 0;
      % inititialize struct cell array to be used to store unphased Q data...
      if isempty(uc.Qpolydata), uc.Qpolydata = cell(1, numel(uc.bas)); end
      for i=1:numel(uc.bas)           % loop over basis set objects in unit cell
        opts.data = uc.Qpolydata{i};  % its allotted poly storage struct
        % call generic discrep routine (which is overloaded for QP LP bases):
        [Qb uc.Qpolydata{i}] = uc.bas{i}.evalunitcelldiscrep(uc, opts);
        Q = [Q Qb];    % stack as cols in basis' order (also works for 3d poly)
      end
    end % func

    function [B B1 B2 d] = evalbasestargetcopies(uc, p, opts) % ..B row from pts
    % EVALBASESTARGETCOPIES - matrix mapping UC bases to val etc at a pointset
    %
    %  Br = EVALBASESTARGETCOPIES(uc, p, opts) evaluates values
    %   opts.nei = 0,1: controls of 1x1 or 3x3 summation over pointset copies
    %   opts.poly: if true, returns 3d array (:,:,1:5) of alpha polynomial mats
    %
    % See also basis/EVALTARGETCOPIES
      if nargin<3, opts = []; end
      if isfield(opts,'data') & numel(opts.data)==numel(uc.bas)
        d = opts.data;
      else
        d = cell(1, numel(uc.bas));   % init storage cell array for this pts
      end
      na = max(1,nargout-1);          % how many matrix output args wanted?
      o = opts;                       % pass through stuff via options
      o.dom = uc;                     % will always be evaluating in unit cell
      [N noff] = uc.setupbasisdofs;
      M = numel(p.x);
      wantpoly = 0; if isfield(opts, 'poly') & opts.poly, wantpoly = 1; end
      if wantpoly, B = zeros(M,N,5); else B = zeros(M,N); end % preallocate
      if na>1, B1 = B; end                                    % (much faster!)
      if na>2, B2 = B; end
      for i=1:numel(uc.bas)           % loop over basis set objects in unit cell
        b = uc.bas{i}; ns = noff(i)+(1:b.Nf);
        o.data = d{i};                % its allotted poly storage struct
        if na==1;                     % get same # args out as requested
          if wantpoly, [B(:,ns,:) d{i}] = b.evaltargetcopies(p, uc, o); 
          else [B(:,ns) d{i}] = b.evaltargetcopies(p, uc, o); end
        elseif na==2
          if wantpoly, [B(:,ns,:) B1(:,ns,:) d{i}] = b.evaltargetcopies(p, uc, o);
          else [B(:,ns) B1(:,ns) d{i}] = b.evaltargetcopies(p, uc, o); end
        elseif na==3
          if wantpoly
            [B(:,ns,:) B1(:,ns,:) B2(:,ns,:) d{i}] = b.evaltargetcopies(p,uc,o);
          else
            [B(:,ns) B1(:,ns) B2(:,ns) d{i}] = b.evaltargetcopies(p, uc, o);
          end
        end
      end
      if na==1, B1 = d; elseif na==2, B2 = d; end % get output arg list correct 
    end % func
    
    function addqpuclayerpots(uc, k, opts) % ...................................
    % ADDQPUCLAYERPOT - adds a sticking-out layer potential basis to unit cell
    %
    %  makes 4 instances of qpuclayerpot and makes them influence the unit cell
      if nargin==1 | isempty(k), k = uc.k; end
      uc.bas = {uc.bas{:}, qpuclayerpot(uc, 'L', 'd', k), qpuclayerpot(uc, 'L', 's', k), qpuclayerpot(uc, 'B', 'd', k), qpuclayerpot(uc, 'B', 's', k)};
    end
    
  end % methods
end
