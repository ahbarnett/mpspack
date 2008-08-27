classdef basis < handle
    
  % Class basis - Abstract class that defines the interfaces which are
  % common for all basis objects. Also defines generic qpunitcell
  % discrepancy evaluation routine.
    
  properties
    k                       % Wavevector
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
  
  methods % ....................................... actual methods
  
    function updateNf(b) % ....................... dummy for all but LP bases
    end
  
    function [Q d] = evalunitcelldiscrep(b, uc, opts)
    % EVALUNITCELLDISCREP - matrix which maps basis coeffs to discrepancy on UC
    %
    %  Default routine for bases which don't know anything special about the UC.
    %
    %  Q = EVALUNITCELLDISCREP(bas, uc, opts) uses options in opts including:
    %   opts.dom: the domain in which the basis is evaluated (eg, conn. dom)
    %   opts.nei = 0,1: 0 has no source nei copies, 1 has 3x3 (used to fill C).
    %   opts.data: if present, treats as struct containing ab-poly matrix data
    %   opts.poly: if true, returns 3d-array Q(:,:,1:5) of polynomial coeff
    %    matrices for alpha^{-2}, alpha^{-1}, ... , alpha^2.
    %
    %  [Q data] = EVALUNITCELLDISCREP(bas, uc, opts) passes out data struct
    %   containing all matrix coeffs needed to recompute Q at a new (alpha,beta)
    %   This struct may then be passed in as opts.data for future fast evals.
    %
    %  To Do: fix nei=1 to return Q-alpha-poly, for filling C !
    %
    % See also QPUNITCELL.
      if ~isa(uc, 'qpunitcell'), error('uc must be a qpunitcell object!'); end
      if nargin<3, opts = []; end
      if isfield(opts, 'nei'), nei = opts.nei; else nei = 0; end
      if nei~=0 & nei~=1, error('opts.nei must be 0 or 1!'); end
      if ~isfield(opts, 'dom'), opts.dom = uc; end     % eval in uc by default
      if isfield(opts, 'data'), d = opts.data; else d = []; end
      recompute = isempty(d);
      if ~recompute              % make sensible guess as to if needs bas evals
        recompute = (opts.dom.k~=d.k | nei~=d.nei | numel(uc.L.x)~=size(d.L,1));
      end
      if recompute, d.k = uc.k; d.nei = nei; end
      wantpoly = 0; if isfield(opts, 'poly') & opts.poly, wantpoly = 1; end
      if nei==0                           % no neighboring local copies
        if recompute
          [d.L d.Ln] = b.eval(uc.L, opts);   % do expensive basis evaluations
          [d.R d.Rn] = b.eval(uc.R, opts);
          [d.B d.Bn] = b.eval(uc.B, opts);
          [d.T d.Tn] = b.eval(uc.T, opts);
        end
        if wantpoly
          ML = size(d.L,1); MB = size(d.B,1); % in case L,R diff # pts from B,T
          Q = zeros(2*(ML+MB), size(d.L,2), 5);
          Q(1:2*ML,:,2) = [-d.R; -d.Rn];           % alpha^{-1} part
          Q(:,:,3) = [d.L; d.Ln; d.B - (1/uc.b)*d.T; d.Bn - (1/uc.b)*d.Tn];
        else
          Q = [d.L - (1/uc.a)*d.R; d.Ln - (1/uc.a)*d.Rn; ...
               d.B - (1/uc.b)*d.T; d.Bn - (1/uc.b)*d.Tn]; % formula for discrep
        end
        
      else           % ******** to fix to include wantpoly=1 (the translates could be replaced with pointset copies translated by hand, for speed too): ************                     % sum local neighbors w/ phases
        nL = 2*numel(uc.L.x); nB = 2*numel(uc.B.x); % dofs in f,f', then g,g'
        Q = zeros(nL+nB, b.Nf);                     % zero the Q matrix
        c = uc.L.translate(2*nei*uc.e1 -(nei+1)*uc.e2); % init loc of seg copy
        for i=-nei:nei
          c.translate(-(2*nei+1)*uc.e1 + uc.e2); % move to left side & up one
          [C Cn] = b.eval(c, opts);
          Q(1:nL,:) = Q(1:nL,:) + (uc.a^nei)*(uc.b^(-i)) * [C; Cn];
          c.translate((2*nei+1)*uc.e1); % move to right side, same row
          [C Cn] = b.eval(c, opts);
          Q(1:nL,:) = Q(1:nL,:) - (uc.a^(-1-nei))*(uc.b^(-i)) * [C; Cn];
        end
        c = uc.B.translate(2*nei*uc.e2 -(nei+1)*uc.e1); % init loc of seg copy
        for i=-nei:nei
          c.translate(-(2*nei+1)*uc.e2 + uc.e1); % do as above, flipped e1-e2
          [C Cn] = b.eval(c, opts);
          Q(nL+1:end,:) = Q(nL+1:end,:) + (uc.b^nei)*(uc.a^(-i)) * [C; Cn];
          c.translate((2*nei+1)*uc.e2);
          [C Cn] = b.eval(c, opts);
          Q(nL+1:end,:) = Q(nL+1:end,:) - (uc.b^(-1-nei))*(uc.a^(-i)) * [C; Cn];
        end
      end
    end % func
    
    function [B B1 B2 d] = evaltargetcopies(b, p, uc, opts) % .. sum over target
    % EVALTARGETCOPIES - eval basis on (maybe sum neighbor copies of) pointset
    %
    %  B = evaltargetcopies(bas, pts, uc, opts)
    %  [B data] = evaltargetcopies(bas, pts, uc, opts)
    %  [B Bn data] = evaltargetcopies(bas, pts, uc, opts)
    %  [B Bx By data] = evaltargetcopies(bas, pts, uc, opts)
    %
    %  opts.poly true makes B (B1, B2 if requested) be alpha-poly matrix coeffs
    %   of size MxNx5, as opposed to the usual MxN.
    %  opts.nei=0 this routine just calls eval. (Apart from it may return matrix
    %   poly, and can store the matrix for later use)
    %  opts.nei=1 it computes the weighted mismatch of eval (using lattice
    %   vectors from uc).
    %  opts.copylist overrides opts.nei and uses an arbitrary set of target copy
    %   translations (copylist.t), alpha-powers (copylist.apow) and remaining
    %   phase factors (copylist.remph), eg see qpuclayerpot/evaltargetcopies
    %
    %  opts.dom gets passed through to basis.eval
    %
    % See also EVAL.
      if ~isa(uc, 'qpunitcell'), error('uc must be a qpunitcell object!'); end
      if nargin<4, opts = []; end
      if isfield(opts, 'nei'), nei = opts.nei; else nei = 0; end
      if nei~=0 & nei~=1, error('opts.nei must be 0 or 1!'); end
      if isfield(opts, 'data'), d = opts.data; else d = []; end       % d=data
      wantpoly = 0; if isfield(opts, 'poly') & opts.poly, wantpoly = 1; end
      na = max(1,nargout-1);          % how many args out of bas.eval 
      N = b.Nf;
      M = numel(p.x);
      if isfield(opts, 'copylist')
        c = opts.copylist;
      else                       % build copy info list for default mismatch
        i = 1;
        for n=-nei:nei, for m=-nei:nei              % loop over 1x1 or 3x3 grid
            c.t(i) = n*uc.e1 + m*uc.e2;             % translations
            c.apow(i) = -n;                         % alpha-powers, and ...
            c.remph(i) = (-1)^n * (-1/uc.b)^m;      % remaining phase fac (beta)
            i = i+1;
          end, end
      end
      nc = numel(c.t);
      recompute = isempty(d);
      %if ~recompute % make sensible guess as to if needs bas evals (slow test)
      %  recompute = (b.k~=d.k | nc~=size(d.B,3) | numel(p.x)~=size(d.B,1) | p~=d.p | (na>1 & ~isfield(d,'B1')) | (na>2 & ~isfield(d,'B2')));
      %end
      if recompute
        d.k = b.k; % records stored k for d.B etc; uses k from basis not uc
        d.p = p;   % record identity of pointset used (recomputes if changes)
        d.B = zeros(M,N,nc); if na>1, d.B1 = zeros(M,N,nc); end  % alloc storage
        if na>2, d.B2 = zeros(M,N,nc); end
        for i=1:nc % could write this loop inside each na choice, faster...
          if c.t(i)==0
            pt = p;                             % preserve identity (jump rels)
          else
            pt = pointset(p.x + c.t(i), p.nx);  % make moved target copy
          end
          if na==1                          % call eval w/ appropriate # args
            d.B(:,:,i) = b.eval(pt, opts);
          elseif na==2
            [d.B(:,:,i) d.B1(:,:,i)] = b.eval(pt, opts);
          else
            [d.B(:,:,i) d.B1(:,:,i) d.B2(:,:,i)] = b.eval(pt, opts);
          end
        end
      end
      if wantpoly                       % ------ alpha matrix polynomials
        B = zeros(size(d.B,1), size(d.B,2), 5);
        for i=1:nc
          B(:,:,3+c.apow(i)) = B(:,:,3+c.apow(i)) + c.remph(i)*d.B(:,:,i);
        end
        if na>1, B1 = zeros(size(B));
          for i=1:nc
            B1(:,:,3+c.apow(i)) = B1(:,:,3+c.apow(i)) + c.remph(i)*d.B1(:,:,i);
          end
        end
        if na>2, B2 = zeros(size(B));
          for i=1:nc
            B2(:,:,3+c.apow(i)) = B2(:,:,3+c.apow(i)) + c.remph(i)*d.B2(:,:,i);
          end
        end
      else                              % ------ just B block wanted, no poly
        ph = uc.a.^c.apow .* c.remph;        % list of phases
        B = ph(1) * d.B(:,:,1);              % sum basis eval mats with phases
        for i=2:nc, B = B + ph(i) * d.B(:,:,i); end
        if na>1, B1 = ph(1) * d.B1(:,:,1);
          for i=2:nc, B1 = B1 + ph(i) * d.B1(:,:,i); end
        end
        if na>2, B2 = ph(2) * d.B2(:,:,1);
          for i=2:nc, B2 = B2 + ph(i) * d.B2(:,:,i); end
        end
      end
      % get output args correct when not 1 or all 4 wanted...
      if na==1, B1 = d; elseif na==2, B2 = d; end
    end %func
    
    function varargout = evalsourcecopiestargetcopies(b, p, uc, opts)
    %
    %
    %  opts.sourcenei = 0,1: size of source copies neighborhood (1x1 or 3x3)
    %  opts.targetnei = 0,1: size of target copies neighborhood (1x1 or 3x3)
    %  opts.dom chooses domain handle that basis is to evaluate in.
      if ~isa(uc, 'qpunitcell'), error('uc must be a qpunitcell object!'); end
      if nargin<4, opts = []; end
      if isfield(opts, 'sourcenei'), snei = opts.sourcenei; else snei = 0; end
      if snei~=0 & snei~=1, error('opts.sourcenei must be 0 or 1!'); end
      if isfield(opts, 'targetnei'), tnei = opts.targetnei; else tnei = 0; end
      if tnei~=0 & tnei~=1, error('opts.targetnei must be 0 or 1!'); end
      if isfield(opts, 'data'), d = opts.data; else d = []; end       % d=data
      wantpoly = 0; if isfield(opts, 'poly') & opts.poly, wantpoly = 1; end
      na = max(1,nargout-1);          % how many args out of bas.eval
      N = b.Nf;
      M = numel(p.x);
      if snei==0
        opts.nei = tnei;
        [varargout{1:nargout}] = b.evaltargetcopies(p, uc, opts);
      else                            % sourcenei = 1
        if tnei==0
          %set up a src set and pass to evaltargetcopies
          
          
        else
          % set up double-size targ set, ditto.
          
          
        end
      end
    end
    
  end % methods
end
