    % QPUCLAYERPOT - setup a sticking-out QP layer potential basis on unit cell
    %
    %  bas = qpuclayerpot(uc, direc, a) creates a set of 6 or 10 copies of
    %   layer densities on either L or B segment of the unit cell object uc,
    %   depending on if direc is
    %   'L' or 'B'. The argument a controls mix of SLP and DLP, as in layerpot.
    %   uc.buffer = 0:  1xQPLP scheme on UC with 6 copies
    %   uc.buffer = 1:  3xQPLP scheme on UC with 10 copies
    %
    %  Note that in fact the LPs are not constructed until evaluation time, but
    %   the setup of locations and phases is done here. Similarly, the phases
    %   used for the copies depend on (alpha, beta) for the unit cell at the
    %   time of evaluation.
    %  Barnett, 8/11/08, entirely rewritten to use copylists 9/8/08, k-free
    %   interface 7/28/09
    %
    % See also: LAYERPOT, QPUNITCELL

classdef qpuclayerpot < handle & basis
  properties
    uc                              % unit cell that it's attached to
    direc                           % 'L' or 'B'
    pseg                            % principal segment ('p', parallel), L or B
    oseg                            % other ('o') segment, B or L
    lp                              % layer-potential object on principal seg
    copylist                        % cell array of target-transl & phase lists 
    a                               % 1-by-2, mixture weights of SLP and DLP
    quad                            % quadrature option (opts.quad in S,D,T)
  end
  
  methods
    function b = qpuclayerpot(uc, direc, a, opts) % ............ constructor
      if nargin<4, opts = []; end
      b.uc = uc;
      % set up copylist displacements (via a,b-powers) for standard bas eval...
      if uc.buffer==0                   % 1xQPLP scheme: 6 segments per basis
        c{1}.apow = [0 0 0 1 1 1];
        c{1}.bpow = [-1 0 1 -1 0 1];
      elseif uc.buffer==1               % 3xQPLP scheme: 10 segments per basis
        % order of 3 cells is tau_{-1}, tau_0, tau_1 (for the DLP in L direc)
        c{1}.apow = [-1 -1 2 2]; c{1}.bpow = [-1 2 -1 2];
        c{2}.apow = [-1 2];      c{2}.bpow = [0 0];
        c{3}.apow = [-1 -1 2 2]; c{3}.bpow = [-2 1 -2 1];
      end
      if direc=='L' | direc=='l'                     % 'L' direction (left)
        b.pseg = uc.L; b.oseg = uc.B; b.direc = 'L';
      elseif direc=='B' | direc=='b'                 % 'B' direction (bottom)
        for j=1:numel(c)                             % swap a<->b powers
          temp = c{j}.apow; c{j}.apow = c{j}.bpow; c{j}.bpow = temp;
        end
        b.pseg = uc.B; b.oseg = uc.L; b.direc = 'B';
      else
        error('invalid direc!');
      end
      for j=1:numel(c)    % set up (target) translations to match a,b-powers...
        c{j}.t = -(uc.e1*c{j}.apow + uc.e2*c{j}.bpow); % key: e1=alpha, e2=beta
        c{j}.remph = ones(size(c{j}.apow));        % phases all 1
      end
      b.copylist = c;
      if ~isnumeric(a)
        switch a
         case {'single', 'S', 's', 'SLP'}
          a = [1 0];
         case {'double', 'D', 'd', 'DLP'}
          a = [0 1];
        end
      end
      b.lp = layerpot(b.pseg, a, opts);  % create a LP basis on p segment
      b.lp.doms = uc;    % make the lp affect the UC (so lp knows k)
      b.a = a;
    end % func
    
    function Nf = Nf(b, opts)  %...............Nf method reads segment # quadr
      Nf = (1+2*b.uc.buffer) * numel(b.pseg.x);
      % # quadr points on 'parallel' segment, 3x if 3x3 QPLP scheme
    end
    
    function showgeom(b, varargin) % ............. crude show discr pts of segs
      t =  b.copylist{:}.t;     % NB get t from all cell arrays
      for j=1:numel(t)
        seg(j) = b.pseg.translate(-t(j));   % make temp segment copies for plot
      end
      domain.showsegments(seg, varargin{:});  % pass opts if present onwards
    end
    
    function varargout = evalunitcellcopies(b, p, ucdummy, opts) % ...overloads
    % EVALUNITCELLCOPIES - basis eval of QP UC layer potentials, w/ data output
    %
    %  See qpuclayerpots/EVALUNITCELLCOPIES for interface.
    %
    %  Refers to current alpha, beta and segments in unit cell, but all other
    %  properties such as translations were frozen in at construction time.
    %  This a wrapper around the evalunitcellcopies method, for qpuclayerpot.
    %  Note that target-movement is used to achieve source move.
    %  The 3rd input argument is ignored, since the basis already contains the
    %  handle of unit cell.
    %
    %  See also qpuclayerpots/EVAL
      if nargin<4, opts = []; end
      if ~isfield(opts, 'poly'), opts.poly = 0; end
      N = b.Nf;                             % cut for speed? (now a method)
      if b.uc.buffer==0       % 1x1 UC case
        opts.copylist = b.copylist{1};   % pass QPUC's copylist to lp summing...
        [varargout{1:nargout}] = b.lp.evalunitcellcopies(p, b.uc, opts);
      else         % buffer>0, not working yet:
        M = numel(p.x);                           % height of each block column
        Ns = N/3;                                 % basis size for each 3 LPs
        B = zeros(M,N,opts.poly+1);               % preallocate
        if nargout>2, B1 = zeros(M,N,opts.poly+1); end
        if nargout>3, B2 = zeros(M,N,opts.poly+1); end
        for j=1:3
          ns = (j-1)*Ns + (1:Ns);
          opts.copylist = b.copylist{j};
          if nargout==1
            B(:,ns,:) = b.lp.evalunitcellcopies(p, b.uc, opts);
          elseif nargout==2
            [B(:,ns,:) d] = b.lp.evalunitcellcopies(p, b.uc, opts);
          elseif nargout==3
            [B(:,ns,:) B1(:,ns,:) d] = b.lp.evalunitcellcopies(p, b.uc, opts);
          else
            [B(:,ns,:) B1(:,ns,:) B2(:,ns,:) d] = b.lp.evalunitcellcopies(p,...
                                                              b.uc, opts);
          end
        end
        varargout(1) = {B};                      % set up output argument list
        if nargout==2, varargout(2) = {d};
        elseif nargout==3, varargout(2:3) = {B1 d};
        else varargout(2:4) = {B1 B2 d}; end
      end
    end % func
    
    function varargout = eval(b, p, varargin) %....std eval (witout data out)
    % EVAL - basis eval of quasi-periodizing unit-cell layer potential copies
    %
    %  See basis/EVAL for interface.
    %
    %  Refers to current alpha, beta and segments in unit cell, but all other
    %  properties such as translations were frozen in at construction time.
    %  This a wrapper around the evalunitcellcopies method, discarding output
    %  data to match eval's output argument format. If there's input opts.data
    %  it will be used for fast rephasing.
    %
    %  See also basis/EVAL, qpuclayerpots/EVALUNITCELLCOPIES
      [varargout{1:nargout} datadummy] = b.evalunitcellcopies(p, b.uc, ...
                                                        varargin{:});
    end % note how the datadummy asks for one more output arg than passed out.
    
    function [Q d] = evalunitcelldiscrep(b, uc, opts) % ....overloads basis
    % EVALUNITCELLDISCREP - return Q matrix mapping a QPUC basis to UC discrep
    %
    %  Same usage as basis/EVALUNITCELLDISCREP, which it overloads.
    %
    %  Includes source quadr wei, ie dofs are density values, as always for LP.
    %  Currently uses no symmetry relations. Stores Bloch-indep interaction
    %   matrices using evalunitcellcopies data cells in d.
    %  Uses unit cell in 2nd argument.
    %  Ignores opts.nei since does not apply to QPUC lp's.
    %  uc.buffer = 0 gives 1x1 scheme, works
    %  uc.buffer = 1 gives 3x3 scheme, which never worked correctly.
    %
    %  Uses only distant interactions + Id, so has spectral convergence.
    %  See testqpuclayerpot.m for convergence plot
    %
    % See also QPUCLAYERPOT, LAYERPOT, QPUNITCELL, basis/EVALUNITCELLDISCREP.
      if nargin<3, opts = []; end
      if ~isfield(opts, 'dom'), opts.dom = uc; end     % eval in uc by default
      if ~isfield(opts, 'poly'), opts.poly = 0; end    % default no poly
      recompute = ~isfield(opts, 'data');
      if ~recompute    % data exists, could be out of date in a variety of ways
        d = opts.data;
        if isempty(d) || size(d{1}.B,1)~=b.Nf || d{1}.k~=opts.dom.k  % NB relop
          recompute=1; end
      end
      if recompute %--------------use uc.buffer to build copylist then fill data
        d = {};
        if isfield(opts,'data'), opts = rmfield(opts,'data'); end % ignore data
        if uc.buffer==0 %............1x1 unit cell discrep (d has 2 cells)
          c = uc.discrepcopylist(-1:1, 1, b.direc, 0);
          d{1} = b.lp.copiesdata(b.pseg, c.t, 2, opts); d{1}.copylist = c;
          c = uc.discrepcopylist(1, 0:1, b.direc, 3);
          d{2} = b.lp.copiesdata(b.oseg, c.t, 2, opts); d{2}.copylist = c;
        elseif uc.buffer==1 %...........3x3 unit cell discrep (d has 16 cells)
          for i=1:7
            p = i-4;      % parallel seg targets: indices of 7 rows (-3:3) 
            c = uc.discrepcopylist(p, 3, b.direc, 0);
            d{i} = b.lp.copiesdata(b.pseg, c.t, 2, opts); d{i}.copylist = c;
          end
          dc = 8;         % other (non-parallel) seg targets, are cells 8-16
          for i=1:3       % loop over eg g_{-1}, g_0, g_1 block rows
            for j=1:3     % loop over eg tau_{-1}, tau_0, tau_1 block cols
              hippow = mod(j+1,3)+1; % [3 1 2], highest parallel power vs j
              c = uc.discrepcopylist(hippow, (4-i)+[0 -3], b.direc, 1+j);
              d{dc} = b.lp.copiesdata(b.oseg, c.t, 2, opts); d{dc}.copylist = c;
              dc = dc+1;
            end
          end
        end
      end
      % ------------------now use the data in d{:} to fill subblocks of Q matrix
      ind3 = ceil(opts.poly/2)+1;  % index to write Id to, match b/evaluccopies
      MP = numel(b.pseg.x); MO = numel(b.oseg.x); % in case L,B neq # quadr pts
      
      if uc.buffer==0 %............1x1 unit cell discrep (d has 2 cells)
        opts.data = d{1};         % data (including copylist) for parallel seg
        [Qp Qpn dummy] = ...
            b.lp.evalunitcellcopies(b.pseg, uc, opts); % get matrices for f, f'
        % discrepancy jump relations due to local layerpot (NB no factor 1/2)...
        Qp(:,:,ind3)  = Qp(:,:,ind3)  + b.a(2)*eye(b.Nf); % (alpha power = 0)
        Qpn(:,:,ind3) = Qpn(:,:,ind3) - b.a(1)*eye(b.Nf); % TODO: use diagind!
        opts.data = d{2};         % data (including copylist) for other targ seg
        [Qo Qon dummy] = ...
            b.lp.evalunitcellcopies(b.oseg, uc, opts); % get matrices for f, f'
      
      elseif uc.buffer==1 %.......3x3 unit cell discrep (d has 16 cells),broken
        % contribs to same-direction target (parallel 'p') seg: ======
        Qp = zeros(3*MP, 3*MP, opts.poly+1); Qpn = Qp;
        % r gives list of row indices stored in d{1:7} needed for each subblock
        r = {[-3 0] -1 [-2 1]; [-2 1] 0 [-1 2]; [-1 2] 1 [0 3]};
        for i=1:3       % loop over eg f_{-1}, f_0, f_1 block rows in Qp (& Qpn)
          msi = (i-1)*MP + (1:MP); % row indices for subblock
          for j=1:3     % loop over eg tau_{-1}, tau_0, tau_1 block cols
            nsi = (j-1)*MP + (1:MP); % col indices for subblock
            opts.data = d{r{i,j}(1)+4};   % get row data indexed by r{i,j}
            [Qp(msi,nsi,:) Qpn(msi,nsi,:) dummy] = ...
                b.lp.evalunitcellcopies(b.pseg, uc, opts);  % write subblock
            if numel(r{i,j})>1          % if 2nd row contribs to subblock
              opts.data = d{r{i,j}(2)+4};   % get row data indexed by r{i,j}
              [Rij Rijn dummy] = b.lp.evalunitcellcopies(b.pseg, uc, opts);
              Qp(msi,nsi,:) = Qp(msi,nsi,:) + Rij;  % add in row's contrib
              Qpn(msi,nsi,:)= Qpn(msi,nsi,:) + Rijn;
            end
          end
        end
        % discrepancy jump relations due to local layerpot (NB no factor 1/2)...
        Qp(:,:,ind3)  = Qp(:,:,ind3)  + b.a(2)*eye(3*MP); % TODO: use diagind!
        Qpn(:,:,ind3) = Qpn(:,:,ind3) - b.a(1)*eye(3*MP);
        % contribs to other-direction target ('o') seg: ============
        Qo = zeros(3*MO, 3*MO, opts.poly+1); Qon = Qo;
        dc = 8;         % start after the parallel data cells; other are 8-16
        for i=1:3       % loop over eg g_{-1}, g_0, g_1 block rows in Qo (& Qon)
          msi = (i-1)*MO + (1:MO); % row indices for subblock
          for j=1:3     % loop over eg tau_{-1}, tau_0, tau_1 block cols
            nsi = (j-1)*MO + (1:MO); % col indices for subblock
            opts.data = d{dc};   % just dump rephased data from cell dc
            [Qo(msi,nsi,:) Qon(msi,nsi,:) dummy] = ...
                b.lp.evalunitcellcopies(b.oseg, uc, opts);  % write subblock
            dc = dc+1;
          end
        end
      end
      if b.direc=='L' % insure f & g are correct way round. TODO: preallocate Q!
        Q = [Qp; Qpn; Qo; Qon];          % f,f' hit by parallel contribs
      else
        Q = [Qo; Qon; Qp; Qpn];          % or, g,g' hit by parallel contribs
      end
    end % func
      
    function showQdatasegs(b, d) % .................. debug show copylists
    % SHOWQDATASEGS - plot all target segs in copylists in basis' data cells
    %
      np = 1 + 6*b.uc.buffer;  % how many parallel data cells
      no = 1 + 8*b.uc.buffer;  % how many other data cells
      if numel(d)~=np+no, error('unexpected number data cells!'); end
      nsegp = 0; nsego = 0;
      for i=1:np               % parallel cells
        transl = d{i}.copylist.t;
        nsegp = nsegp + numel(transl);
        for j=1:numel(transl)
          b.pseg.translate(transl(j)).plot; hold on;
        end
      end
      for i=np+(1:no)          % other cells
        transl = d{i}.copylist.t;
        nsego = nsego + numel(transl);
        for j=1:numel(transl)
          b.oseg.translate(transl(j)).plot; hold on;
        end
      end
      fprintf('np=%d no=%d, nsegp=%d nsego=%d\n', np, no, nsegp, nsego)
    end % func
      
  end % methods
end
