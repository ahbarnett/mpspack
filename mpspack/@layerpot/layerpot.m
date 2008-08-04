classdef layerpot < handle & basis

% LAYERPOT - create a layer potential basis set on a segment
%
%  b = LAYERPOT(seg, a, k, opts) where a = 'single' or 'double' creates a layer
%   potential basis object with density on segment seg, with wavenumber k
%   (which may be [] in which k is not defined), and options:
%   opts.real: if true, real valued (Y_0), otherwise complex (H_0 outgoing).
%   If a is instead a 1-by-2 array, then a mixture of a(1) single plus a(2)
%   double is created, useful for Brakhage, Werner, Leis & Panich type
%   representations.
%
%   Note that the discretization of the layerpot is given by that of the seg,
%    apart from periodic segments where new quadrature weights may be used.
%
%  To do: real=1 case?

  properties
    real                            % true if fund sol is Y_0, false for H_0^1
    seg                             % handle of segment on which density sits
    a                               % 1-by-2, mixture weights of SLP and DLP
  end

  methods
    function b = layerpot(seg, a, k, opts) %........................ constructor
      if nargin<4, opts = []; end
      if nargin>2 & ~isempty(k), b.k = k; end
      b.N = numel(seg.x);
      b.Nf = b.N;
      if ~isfield(opts, 'real'), opts.real = 0; end
      b.real = opts.real;
      if ~isa(seg, 'segment'), error('seg must be a segment object!'); end
      b.seg = seg;
      if ~isnumeric(a)
        switch a
         case {'single', 'S', 's', 'SLP'}
          a = [1 0];
         case {'double', 'D', 'd', 'DLP'}
          a = [0 1];
        end
      end
      if size(a)~=[1 2], error('a argument is not size 1-by-2!'); end
      b.a = a;
    end % func
    
    function [A Ax Ay] = eval(b, p, o) % .........basis evaluator at points p
    % EVAL - evaluate layer potential on pointset or on a segment, with jump rel
    %
    %  A = EVAL(bas, p, opts)           p is the target pointset or segment
    %  [A An] = EVAL(bas, p, opts)
    %  [A Ax Ay] = EVAL(bas, p, opts)
    %
    %   The optional argument opts may contain the following fields:
    %    opts.layerpotside = +1 or -1: determines from which side of a segment
    %     the limit is approached, for jump relation. +1 is the normal side.
      if nargin<3, o = []; end
      self = (p==b.seg); % if true, local eval on segment carrying density
      if self & o.layerpotside~=1 & o.layerpotside~=-1
        error('opts.layerpotside must be +1 or -1 for self-interaction');
      end
      if self, p = []; end              % tell S, D, etc to use self-interaction
      if nargout==1 %------------------------------- values only
        if b.a(2)==0             % only SLP
          A = layerpot.S(b.k, b.seg, p, o);
          if b.a(1)~=1, A = A * b.a(1); end 
        elseif b.a(1)==0         % only DLP
          A = layerpot.D(b.k, b.seg, p, o);
          if b.a(2)~=1, A = A * b.a(2); end 
        else                     % mixture of SLP + DLP
          A = b.a(1) * layerpot.S(b.k, b.seg, p, o) + ...
              b.a(2) * layerpot.D(b.k, b.seg, p, o);
        end
        if self & b.a(2)~=0
          A = A + o.layerpotside * b.a(2) * eye(size(A)) / 2;  % DLP val jump
        end
        
      elseif nargout==2 %------------------------------- values + normal derivs
        if b.a(2)==0             % only SLP
          A = layerpot.S(b.k, b.seg, p, o);
          o.derivSLP = 1;
          Ax = layerpot.D(b.k, b.seg, p, o);
          if b.a(1)~=1, A = A * b.a(1); Ax = Ax * b.a(1); end 
        elseif b.a(1)==0         % only DLP
          A = layerpot.D(b.k, b.seg, p, o);
          Ax = layerpot.T(b.k, b.seg, p, o);
          if b.a(2)~=1, A = A * b.a(2); Ax = Ax * b.a(2); end 
        else                     % mixture of SLP + DLP
          [A o.Sker] = layerpot.S(b.k, b.seg, p, o);
          [AD o.Dker_noang] = layerpot.D(b.k, b.seg, p, o);
          A = b.a(1) * A + b.a(2) * AD;
          clear AD
          o.derivSLP = 1;
          [Ax o.Dker_noang] = layerpot.D(b.k, b.seg, p, o);
          Dn = layerpot.T(b.k, b.seg, p, o);
          Ax = b.a(1) * Ax + b.a(2) * Dn;
        end
        if self & b.a(2)~=0
          A = A + o.layerpotside * b.a(2) * eye(size(A)) / 2;   % DLP val jump
        end
        if self & b.a(1)~=0
          An = An - o.layerpotside * b.a(1) * eye(size(A)) / 2; % SLP deriv jump
        end
        
      else % ------------------------------------- values, x- and y-derivs
        if b.a(2)==0             % only SLP
          A = layerpot.S(b.k, b.seg, p, o);
          o.derivSLP = 1;
          nx_copy = p.nx; p.nx = ones(size(p.nx)); % overwrite pointset normals
          [Ax o.Dker_noang] = layerpot.D(b.k, b.seg, p, o);
          p.nx = 1i*ones(size(p.nx));
          Ay = layerpot.D(b.k, b.seg, p, o); % reuses H1 (Dker) mat
          p.nx = nx_copy;                          % restore pointset normals
          if b.a(1)~=1, A = A * b.a(1); Ax = Ax * b.a(1); Ay = Ay * b.a(1); end 
       elseif b.a(1)==0         % only DLP
          [A o.Dker_noang] = layerpot.D(b.k, b.seg, p, o);
          nx_copy = p.nx; p.nx = ones(size(p.nx)); % overwrite pointset normals
          [Ax o.Sker] = layerpot.T(b.k, b.seg, p, o);
          p.nx = 1i*ones(size(p.nx));
          Ay = layerpot.T(b.k, b.seg, p, o);
          p.nx = nx_copy;                          % restore pointset normals
          if b.a(2)~=1, A = A * b.a(2); Ax = Ax * b.a(2); Ay = Ay * b.a(2); end
        else                    % mixture of SLP + DLP
          [A o.Sker] = layerpot.S(b.k, b.seg, p, o);
          [AD o.Dker_noang] = layerpot.D(b.k, b.seg, p, o);
          A = b.a(1) * A + b.a(2) * AD;
          clear AD
          o.derivSLP = 1;
          nx_copy = p.nx; p.nx = ones(size(p.nx)); % overwrite pointset normals
          Ax = layerpot.D(b.k, b.seg, p, o);
          Ax = b.a(1) * Ax + b.a(2) * layerpot.T(b.k, b.seg, p, o);
          p.nx = 1i*ones(size(p.nx));
          Ay= layerpot.D(b.k, b.seg, p, o);
          Ay = b.a(1) * Ay + b.a(2) * layerpot.T(b.k, b.seg, p, o);
          p.nx = nx_copy;                          % restore pointset normals
        end
        if self
          error('LP self jump rels not yet implemented for Ax Ay grad vec!')
        end
        
      end
    end % func
  end % methods
  
  methods (Static)
    [A Sker] = S(k, s, t, o)                      % Phi kernel matrix
    [A Dker_noang cosker] = D(k, s, t, o)         % dPhi/dny or dPhi/dny matrix
    [A Sker Dker_noang] = T(k, s, t, o)           % d^2Phi/dnxdny matrix
  end % methods
end