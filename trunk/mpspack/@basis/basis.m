classdef basis < handle
    
  % Class basis - Abstract class that defines the interfaces which are
  % common for all basis objects. Also defines qpunitcell
  % basis and discrepancy evaluation methods that apply to generic bases.
    
  properties
    k                       % Wavevector. TODO: make this refer to a domain's k
    N                       % Number or degree of basis fct. Exact
                            % specification depends on type of basis
    Nf                      % Actual number of basis functions
  end
  
  methods (Abstract)
    [A, A1, A2] = eval(b, pts) 
                             % Evaluate a basis on a set of points
                             % A=eval(pts) returns only fct. values
                             % [A A1]=eval(pts) returns fct. values plus normal
                             % derivatives.
                             % [A A1 A2]=eval(pts) returns fct. values plus
                             % x and y derivatives.
  end
  
  methods % ....................................... actual methods ...........
  
    function updateNf(b) % ................ dummy for all but layerpot bases
    end
  
    function d = copiesdata(b, p, transl, nargs, opts) %.......eval copies data
    % COPIESDATA - data struct of basis func evals on many copies of pointset
    %
    % data = COPIESDATA(bas, pts, transl, nargs)
    %   Computes data struct containing raw (unphased) basis evaluation matrices
    %   over copies. From this data, basis.evalunitcellcopies can rapidly
    %   sum and phase the matrices.
    %   transl is list of translations t (as C-#s).
    %   nargs = 1,2,3: controls whether data struct stores vals, or vals+normal
    %    derivs, or vals+x,y-derivs
    %
    % data = COPIESDATA(bas, pts, transl, nargs, opts) allows options to
    %   be passed to basis.eval, such as:
    %    opts.dom: domain for basis evaluation.
    %          
      if nargin<5, opts = []; end
      nc = numel(transl);
      N = b.Nf;
      M = numel(p.x);
      d.k = b.k; % records stored k for d.B etc; uses k from basis not UC
      d.p = p;   % record identity of pointset used (recomputes if changes)
      d.B = zeros(M,N,nc); if nargs>1, d.B1 = zeros(M,N,nc); end % alloc storage
      if nargs>2, d.B2 = zeros(M,N,nc); end
      for i=1:nc % could write this loop inside each na choice, faster...
        if transl(i)==0
          pt = p;                             % preserve identity (jump rels)
        else
          pt = pointset(p.x + transl(i), p.nx); % make moved target copy, as pts
        end
        % following is a hack specific to layerpot basis, shifting origin...
        if isfield(opts, 'Jfilter'), opts.Jfilter.origin = transl(i); end
        if nargs==1                % call eval w/ appropriate # args
          d.B(:,:,i) = b.eval(pt, opts);
        elseif nargs==2
          [d.B(:,:,i) d.B1(:,:,i)] = b.eval(pt, opts);
        else
          [d.B(:,:,i) d.B1(:,:,i) d.B2(:,:,i)] = b.eval(pt, opts);
        end
      end
    end %func
    
    function [B B1 B2 d] = evalunitcellcopies(b, p, uc, opts) % sum over targets
    % EVALUNITCELLCOPIES - sum basis funcs over neighboring copies of pointset
    %
    %  B = evalunitcellcopies(bas, pts, uc, opts) returns in B a MxNx1 3d array
    %   which when squeezed down to MxN is the matrix giving
    %   the sum of each basis function value at each point in pointset pts.
    %   Each column of squeezed B corresponds to one basis func, in basis set
    %   object bas. uc is the unit cell object.
    %   The type of summation is controlled by opts.nei, which gives the size
    %   of a square block of copies (1x1, 3x3, etc), or if opts.copylist is
    %   present this allows explicit control of copy locations in terms of
    %   integer multiples of the unit cell vectors.
    %
    %  [B data] = evalunitcellcopies(bas, pts, uc, opts) also returns a data
    %   structure which if passed back in via opts.data allows rapid
    %   rephasing without recomputing any basis functions.
    %    data has the following fields:
    %     k (wavenumber), p (handle of pointset), B (& optionally B1, B2, data),
    %     copylist (struct describing the locations and phases of matrix data
    %     in B, B1, B2).
    %
    %  [B Bn data] = evalunitcellcopies(bas, pts, uc, opts) also returns normal
    %   derivatives of the basis functions.
    %
    %  [B Bx By data] = evalunitcellcopies(bas, pts, uc, opts) instead returns
    %   x- and y-derivatives of the basis functions.
    %
    %  opts.nei=0 this routine just calls eval. (Apart from it may return matrix
    %   poly, and can store the matrix for later use)
    %  opts.nei=1,2,.. sums over 3x3, or 5x5, etc, sources, using lattice
    %   vectors and phases from UC. This is done using a target copylist.
    %  opts.copylist overrides opts.nei and uses an arbitrary set of target copy
    %   translations (copylist.t), alpha,beta-powers (copylist.apow, bpow)
    %   and remaining phase factors (copylist.remph).
    %  opts.poly>0 makes B (B1, B2 if requested) be alpha-poly matrix coeffs
    %   of size MxNx(poly+1), as opposed to the usual MxNx1. See code for how
    %   the alpha-powers are determined (there is an offset in indexing).
    %  opts.data: if present, overrides opts.nei and opts.copylist, and reuses
    %   precomputed basis eval copies data allowing rapid evaluation. 
    %
    %  opts.dom gets passed through to basis.eval
    %
    % See also EVAL, COPIESDATA
      if ~isa(uc, 'qpunitcell'), error('uc must be a qpunitcell object!'); end
      if nargin<4, opts = []; end
      if ~isfield(opts, 'poly'), opts.poly = 0; end       % default no poly
      na = max(1,nargout-1);            % how many args wanted out of bas.eval?
      recompute = ~isfield(opts, 'data');
      % insert more detailed test for if data is current... (quick for now)
      %if ~recompute % make sensible guess as to if needs bas evals (slow test)
      %  recompute = (b.k~=d.k | nc~=size(d.B,3) | numel(p.x)~=size(d.B,1) | p~=d.p | (na>1 & ~isfield(d,'B1')) | (na>2 & ~isfield(d,'B2'))); % TODO: speed up!
      if ~recompute
        d = opts.data; if isempty(opts.data) | (b.k~=d.k), recompute=1; end
      end
      if recompute
        d = [];
        if ~isfield(opts, 'copylist')    % make default square local copy grid
          if isfield(opts, 'nei'), nei = opts.nei; else nei = 0; end
          i = 1;
          for n=-nei:nei, for m=-nei:nei            % loop over 3x3, etc, grid
              c.t(i) = -n*uc.e1 - m*uc.e2;          % translate (-ve of src)
              c.apow(i) = n; c.bpow(i) = m;         % alpha,beta-powers, and ...
              c.remph(i) = 1;                       % any remaining phase fac
              i = i+1;
            end, end
        else
          c = opts.copylist;                        % use passed-in copylist
        end
        d = copiesdata(b, p, c.t, na, opts); % Do the intensive basis evals!
        d.copylist = c;                      % since copiesdata doesn't write
      else
        c = d.copylist;                      % either way, c now has copylist
      end
      B = zeros(size(d.B,1), size(d.B,2), opts.poly+1); % alloc poly size + 1
      if opts.poly>0                        % want poly coeffs
        ind3 = ceil(opts.poly/2)+1+c.apow;  % third index to write to (NB shift)
        if max(ind3)>opts.poly+1 | min(ind3)<1  % catch overflow in ind3
          error('Not enough alpha matrix poly powers avail! (try increasing opts.poly)');
        end
        phase = uc.b.^c.bpow.*c.remph;      % phases, don't include a-powers
       else
        ind3 = ones(size(c.apow));          % if no poly, third index always 1
        phase = uc.a.^c.apow.*uc.b.^c.bpow.*c.remph;       % include a,b powers
      end
      nc = numel(c.t);                      % # copy contribs
      for i=1:nc   % ... loop over copy contributions, adding in (indep of poly)
        B(:,:,ind3(i)) = B(:,:,ind3(i)) + phase(i)*d.B(:,:,i);
      end
      if na>1, B1 = zeros(size(B));
        for i=1:nc
          B1(:,:,ind3(i)) = B1(:,:,ind3(i)) + phase(i)*d.B1(:,:,i);
        end
      end
      if na>2, B2 = zeros(size(B));
        for i=1:nc
          B2(:,:,ind3(i)) = B2(:,:,ind3(i)) + phase(i)*d.B2(:,:,i);
        end
      end
      % get output args correct when neither 1 nor all 4 args wanted...
      if na==1, B1 = d; elseif na==2, B2 = d; end
    end %func
    
    function [C d] = evalunitcelldiscrep(b, uc, opts) % ........... discrepancy
    % EVALUNITCELLDISCREP - matrix mapping basis coeffs to discrepancy on (3x)UC
    %
    %  Generic routine for bases which don't know anything special about the UC.
    %
    %  C = EVALUNITCELLDISCREP(bas, uc, opts) uses options in opts including:
    %   opts.dom: the domain in which the basis is evaluated (eg, conn. dom)
    %   opts.nei = 0,1,2: 1x1, 3x3, or 5x5 source copies (used to fill C).
    %   opts.data: if present, treats as cell array of copies data structs
    %   opts.poly: if true, returns 3d-array Q(:,:,:) of alpha-polynomial coeff
    %    matrices with size(Q,3)=poly+1 and powers given by evalunitcellcopies.
    %
    %  [C data] = EVALUNITCELLDISCREP(bas, uc, opts) passes out data struct
    %   containing all matrix coeffs needed to recompute C at a new (alpha,beta)
    %   This struct may then be passed in as opts.data for future fast evals.
    %
    % See also QPUNITCELL.
      if ~isa(uc, 'qpunitcell'), error('uc must be a qpunitcell object!'); end
      if nargin<3, opts = []; end
      if ~isfield(opts, 'poly'), opts.poly = 0; end       % default no poly
      if isfield(opts, 'nei'), nei = opts.nei; else nei = 0; end
      recompute = ~isfield(opts, 'data');
      % TODO:insert more detailed test for if data is current... (quick for now)
      % recompute = (opts.dom.k~=d.k | nei~=d.nei |ucbuf~=d.ucbuf); %| numel(uc.L.x)~=size(d.L,1));
      if ~recompute
        d = opts.data; if isempty(opts.data) | (b.k~=d{1}.k), recompute=1; end
      end
      if recompute  % use nei and uc.buffer to build copylist then fill data...
        d = {};
        if isfield(opts,'data'), opts = rmfield(opts,'data'); end % ignore data
        if ~isfield(opts, 'dom'), opts.dom = uc; end     % eval in uc by default
        dc = 1;                     % counter over d data cells
        for direc='LB' % ------- meas f then g, filling the data array d{:}
          if direc=='L', seg = uc.L; else seg = uc.B; end  % set parallel seg
          if uc.buffer==0 %............1x1 unit cell discrep (d has 2 cells)
            if nei==0                           % 1x1 src, simplest discrep
              c = uc.discrepcopylist(0, 0, direc, 1);
            elseif nei==1                       % 3x3 source copies w/ cancel
              c = uc.discrepcopylist(-1:1, 1, direc, 1);
            else
              error('recompute data: nei>1 not supported for uc.buffer=0');
            end
          elseif uc.buffer==1 %...........3x3 unit cell discrep (d has 6 cells)
            if nei<2                           % use no source cancellation
              
            elseif nei==2                      % the full scheme, cancellation
              
            else
              error('recompute data: only nei=2 supported for uc.buffer=1');
            end
          elseif uc.buffer==2 %...........3x rescaled unit cell walls
            seg = seg.scale(3);                % new scaled, NOT shifted, seg
            if nei<2                           % use no source cancel
              c = uc.discrepcopylist(-nei:nei, -nei:nei, direc, -3);
            elseif nei==2                      % cancellation leaving 3x5 blocks
              c = uc.discrepcopylist(-nei:nei, 0:2, direc, -5);
            else
              error('recompute data: nei>2 not supported for uc.buffer=2');
            end            
          end
          % calc basis evals A,An for seg target copy list, put in d{dc}...
          d{dc} = b.copiesdata(seg, c.t, 2, opts);
          d{dc}.copylist = c;
          dc = dc+1;
        end
      end
      % -------------now use the data in d{:} to fill subblocks of C matrix...
      % Note this code is now agnostic as to whether poly is wanted or not...
      ML = numel(uc.L.x); MB = numel(uc.B.x);   % in case L diff # pts from B
      N = b.Nf;                                 % number of columns in C
      if uc.buffer~=1 %............1x1 unit cell discrep (any nei); 3x scaled
        C = zeros(2*(ML+MB), N, opts.poly+1);
        opts.data = d{1};         % data (including copylist) for L parallel
        [C(1:ML,:,:) C(ML+(1:ML),:,:) dummy] = ...
            b.evalunitcellcopies(uc.L, uc, opts); % matrix for f, f'; uc.L dummy
        opts.data = d{2};         % data (including copylist) for B parallel
        [C(2*ML+(1:MB),:,:) C(2*ML+MB+(1:MB),:,:) dummy] = ...
            b.evalunitcellcopies(uc.B, uc, opts); % matrix for g, g'; uc.B dummy
      else %........................3x3 unit cell discrep full scheme
        
      end
    end % func
    
  end % methods
end
