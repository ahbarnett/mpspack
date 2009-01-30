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
    recip                % 2x2 matrix with columns the reciprocal lattice vecs
    L, B, R, T           % L, B, R, T segments
    N                    % total # degrees of freedom in QP basis sets
    basnoff              % dof offsets for QP basis sets (numeric list)
    Qpolydata            % matrix polynomial storage for each bas (cell array)
    buffer               % 0,1: simple discrepancy, or 3x3 unit cell discrepancy
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
      uc.L = L; uc.B = B;                     % note LBRT normals neq domain's!
      uc.recip = 2*pi*inv([real(e1) real(e2); imag(e1) imag(e2)]'); % cols
      uc.e1 = e1; uc.e2 = e2;
      uc.r1 = uc.recip(1,1) + 1i*uc.recip(2,1);   % reciprocal lattice vecs
      uc.r2 = uc.recip(1,2) + 1i*uc.recip(2,2);
      uc.k = k; uc.setbloch;                      % default Bloch vector
      uc.buffer = 0;                              % default discrep UC size
    end % func
    
    function [x i j] = fold(uc, x)
    % FOLD - fold a point back into the principal unit cell
    %
    %  [x i j] = FOLD(x) returns in x the folded C-number coordinates of a
    %   list of numbers, and (optionally)
    %   in i and j the integer numbers of unit cells
    %   needed to move in.
      M = [real(uc.e1) real(uc.e2); imag(uc.e1) imag(uc.e2)];
      xvecs = [real(x(:))'; imag(x(:))']; % convert to 2-by-n stack of vectors
      y = M\xvecs + 0.5;                  % work relative to lower left corner
      f = floor(y);
      i = reshape(f(1,:), size(x)); j = reshape(f(2,:), size(x));
      y = M*(y - f - 0.5);
      x = reshape(y(1,:) + 1i*y(2,:), size(x));  % convert from vec to C-number
    end
    
    function [N noff] = setupbasisdofs(uc)
    % SETUPBASISDOFS - set up indices of basis degrees of freedom in unit cell
    %
    %  Warning: slow, taking 0.3 ms even if only one basis set!
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
    
    function h = plot(uc, o)   % ..................... overloads domain plot
    % PLOT - paired-down version of DOMAIN/PLOT for unit-cell objects
      if nargin<2, o = []; end
      h = domain.showsegments(uc.seg, 1, o);          % show 4 segments
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
      text(real(uc.r1)/2, imag(uc.r1)/2, 'X', ...
           'fontsize', 12, 'Fontweight', 'demi');
      text(real(uc.r1+uc.r2)/2, imag(uc.r1+uc.r2)/2, ...
           'M', 'fontsize', 12, 'Fontweight', 'demi');
    end % func
    
    function [B B1 B2 d] = evalbaseswithdata(uc, p, opts) % evalbases w/ data IO
    % EVALBASESWITHDATA - UC version of domain/EVALBASES w/ prestored data I/O
    %
    % B = EVALBASESWITHDATA(uc, pts) evaluates matrix mapping basis coeffs in
    %   the unit cell uc to their values on the pointset pts. This is identical
    %   to a generic uc.evalbases(pts) call, so there's not much point in that.
    %
    % B = EVALBASESWITHDATA(uc, pts, opts) allows options including:
    %   opts.data: cell array of data structs, one for each basis in UC, allows
    %    rapid resummation of prestored submatrices for basis eval with phases
    %    given in uc.
    %
    % [B d] = EVALBASESWITHDATA(...) also returns d a cell array of data struct
    %   of the same form as opts.data may have passed in.
    % [B Bn d] = etc
    % [B Bx By d] = etc
    %
    % See also domain/EVALBASES
      if nargin<3, opts = []; end
      if ~isfield(opts, 'poly'), opts.poly = 0; end
      o = opts;          % opts copy to pass through, but with data overwritten
      if ~isfield(opts, 'data') | isempty(opts.data)
        opts.data = cell(1, numel(uc.bas)); end         % initialize empty data
      o.dom=uc; o.nei=0;              % always eval in the uc domain, no nei blk
      M = numel(p.x);                 % height of each block column
      N = uc.N; noff = uc.basnoff;    % since uc.setupbasisdofs takes 0.3ms
      B = zeros(M,N,opts.poly+1);                          % preallocate
      if nargout>2, B1 = zeros(M,N,opts.poly+1); end
      if nargout>3, B2 = zeros(M,N,opts.poly+1); end
      for i=1:numel(uc.bas) %-----------loop over basis set objects in unit cell
        b = uc.bas{i}; ns = noff(i)+(1:b.Nf);  % column indices for this basis
        o.data = opts.data{i};              % its allotted poly storage struct
        if nargout==1                       % handle each # outp args separately
          Bb = b.evalunitcellcopies(p, uc, o);  % block-column of B
          B(:,ns,:) = Bb;
        elseif nargout==2
          [Bb db] = b.evalunitcellcopies(p, uc, o);
          B(:,ns,:) = Bb; d{i} = db;          % create d cell array of structs
        elseif nargout==3
          [Bb B1b db] = b.evalunitcellcopies(p, uc, o);
          B(:,ns,:) = Bb; B1(:,ns,:) = B1b; d{i} = db;
        else
          [Bb B1b B2b db] = b.evalunitcellcopies(p, uc, o);
          B(:,ns,:) = Bb; B1(:,ns,:) = B1b; B2(:,ns,:) = B2b; d{i} = db;
        end
      end                  %------------
      % get output args correct when neither 1 nor all 4 args wanted...
      if nargout==2, B1 = d; elseif nargout==3, B2 = d; end
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
    %   opts.poly: if true, returns 3d-array Q(:,:,poly+1) of polynomial coeff
    %    matrices of order determined in basis.evalunitcellcopies
      if nargin<2, opts = []; end
      if ~isfield(opts, 'poly'), opts.poly = 0; end
      N = uc.N; noff = uc.basnoff;    % since uc.setupbasisdofs takes 0.3ms
      M = (2 + 4*uc.buffer)*(numel(uc.L.x)+numel(uc.B.x)); % # discrep vals
      Q = zeros(M,N,opts.poly+1);                          % preallocate
      % inititialize struct cell array to be used to store unphased Q data...
      if isempty(uc.Qpolydata), uc.Qpolydata = cell(1, numel(uc.bas)); end
      for i=1:numel(uc.bas)           % loop over basis set objects in unit cell
        b = uc.bas{i}; ns = noff(i)+(1:b.Nf);
        opts.data = uc.Qpolydata{i};  % its allotted poly storage struct
        % call generic discrep routine (which is overloaded for QP LP bases):
        [Q(:,ns,:) uc.Qpolydata{i}] = b.evalunitcelldiscrep(uc, opts);
      end
    end % func

    function addqpuclayerpots(uc, k, opts) % ...................................
    % ADDQPUCLAYERPOT - adds a sticking-out layer potential basis to unit cell
    %
    %  makes 4 instances of qpuclayerpot and makes them influence the unit cell
      if nargin==1 | isempty(k), k = uc.k; end
      uc.bas = {uc.bas{:}, qpuclayerpot(uc, 'L', 'd', k), qpuclayerpot(uc, 'L', 's', k), qpuclayerpot(uc, 'B', 'd', k), qpuclayerpot(uc, 'B', 's', k)};
      uc.setupbasisdofs;
    end
    
    function c = discrepcopylist(uc, ppows, opows, direc, refl)
    % DISCREPCOPYLIST - build segment copylist (rect block + mirror image) for C
    %
    %  This is for generic basis objects that do not know anything about the UC.
    %  The lists returned here are product-grid blocks of copies, minus their
    %   (integer) reflections about half-way across the unit cell. Such lists
    %   are useful building blocks to build the C matrix for generic inclusion
    %   basis sets via basis.evalunitcelldiscrep.
    %
    %  Eg, if ppows=[-1 2], opows=[1 3], refl=0, then 8 copies will result with
    %   p-powers=[-1 -1 2 2 -1 -1 2 2] and o-powers=[1 3 1 3 -1 -3 -1 -3].
    %   The last four of these will have a -1 phase (remph), the first four +1.
    %   Finally, if direc='L' then p-powers is beta power, o-powers is alpha.
    % 
    %  copylist = DISCREPCOPYLIST(uc, ppows, opows, direc, refl) where
    %   direc is 'L' or 'B' returns copylist struct as needed for
    %   bas.evalunitcelldiscrep (refl=1), or qpuclayerpots.evalunitcelldiscrep
    %   (refl=0,2). Other inputs are:
    %   ppowrange = integer list of powers in parallel direction.
    %   opowrange = integer list of powers in other direction.
    %   refl controls what duplication of rect copy array is performed:
    %    refl=0: reflection about line of principal segment
    %            (ie negate the 'o' power)
    %    refl=1: reflection about line parallel to principal segment through
    %            UC center (ie negate the 'o' power and subtract 1)
    %    refl=2: reflection about line of `other' segment
    %            (ie negate the 'p' power)
    %    refl=3: reflection about line parallel to `other' segment through
    %            UC center (ie negate the 'p' power and subtract 1)
    %    refl=4: reflection about line parallel to `other' segment through
    %            UC top (ie negate the 'p' power and subtract 2)
    %    refl<0: duplicate of block translated refl in the 'other' direction
    %            (phase is still negated in the duplicate)
      [ppow opow] = meshgrid(ppows, opows);  % make product grid
      ppow = ppow(:); opow = opow(:);
      if refl==0 | refl==1
        opowr = [opow; -refl-opow]; ppowr = [ppow; ppow];
      elseif refl>=2
        opowr = [opow; opow]; ppowr = [ppow; -(refl-2)-ppow];  % cases refl=2-4
      else                                                     % refl<0
        opowr = [opow; opow+refl]; ppowr = [ppow; ppow];
      end
      if direc=='L' | direc=='l'
        c.t = -uc.e1*opowr - uc.e2*ppowr;   % - signs since translate opp power
        c.apow = opowr; c.bpow = ppowr;
      elseif direc=='B' | direc=='b'
        c.t = -uc.e2*opowr - uc.e1*ppowr;
        c.apow = ppowr; c.bpow = opowr;
      else
        error('uc.discrepcopylist: bad direc');
      end
      c.remph = [ones(size(ppow)); -ones(size(ppow))]; % -ve in discrep form
    end % func
    
    function f = spuriousdistance(uc)
    % SPURIOUSDISTANCE - function tries to kill spurious singvals in 3xUC scheme
    % uses k overall wavenumber, and alpha, beta, from uc.
      om = uc.k;                       % omega = overall wavenumber
      n=ceil(3*om/min(svd(uc.recip))); % size of rect array of points in k-space
      [xx yy] = meshgrid((-n:n)/3);    % now kill k's at centers of recip UCs
      j = find(xx~=round(xx) | yy~=round(yy)); ks = xx(j)*uc.r1 + yy(j)*uc.r2;
      %figure; plot(real(ks), imag(ks), '+'); axis equal;
      % we are left with ks a list, punctured grid of centers in C plane
      d = abs(uc.kbloch - ks);         % list of distances of k from centers
      f = min(abs(d-om));              % min dist to omega-circles at centers
    end
    
    d = datawrapR(uc, d, i)             % trial wrapping-over-R-wall of mat data
    
  end % methods
end
