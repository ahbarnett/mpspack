classdef ftylayerpot < handle & basis

% FTYLAYERPOT - create a Fourier-transform layer potential basis in y-direc
%
%  b = FTYLAYERPOT(orig, a, opts) where a = 'single' or 'double' creates a
%   fourier layer potential basis object with density on vertical ray through
%   the point orig. This means a layer potential whose degrees of freedom are
%   at quadrature points in the 1D Fourier variable along the ray direction.
%   As usual for basis sets, the wavenumber k is
%   determined by that of the affected domain.
%   Sommerfeld contour is used in complex k_y plane. 
%
%   Options: opts.M sets # of quadrature points on Sommerfeld contour.
%   Other options: see requadrature.

% Copyright (C) 2008-2010, Alex Barnett
  properties
    orig                            % location of ray origin (x, and y shift)
    quad, om, si, maxy, minx, shift % quad, om, nearsing(=si), etc, Som. contour
    a                               % amounts of single and double layer
    kj, wj                          % k_y quadr nodes and weights
  end

  methods
    function b = ftylayerpot(orig, a, opts) %...................... constructor
      if nargin<1, orig = 0; end               % default (needed so can copy)
      if nargin<2, a = 's'; end                % default
      if nargin<3, opts = []; end
      if ~isfield(opts, 'M'), opts.M = []; end
      if numel(orig)~=1, error('orig must be a single number!'); end
      b.orig = orig;
      if ~isnumeric(a)
        switch a
         case {'single', 'S', 's', 'SLP'}
          a = [1 0];
         case {'double', 'D', 'd', 'DLP'}
          a = [0 1];
         otherwise
          error('cannot interpret non-numeric input argument a!');
        end
      end
      if size(a)~=[1 2], error('a argument is not size 1-by-2!'); end
      b.a = a;
      if ~isfield(opts,'omega'), opts.omega = 1; end    % default freq
      b.nmultiplier = 1;                                % default
      b.requadrature(opts.M, opts);     % construct quadrature points...
    end % func
    
    function Nf = Nf(b)  % ............... Nf method reads # quadr pts
      Nf = numel(b.kj);
    end

    function requadrature(ba, N, opts)
    % REQUADRATURE - make N-pt Sommerfeld contour quadrature for a ftlayerpot.
    %
    % Uses freq omega from the first affected domain.
    % For curve details, see Hankel_integral.m
    %   Options:
    %    opts.omega = overrides freq for choosing Sommerfeld quadrature contour
    %    opts.nearsing = distance to nearest singularity on Re axis.
    %    opts.shift = translation applied to all kj quadr pts
    %    opts.maxy = largest vertical location for evaluation (controls quadr)
    %    opts.minx = smallest x distance for evaluation (controls quadr)
    %    opts.quad = overrides quadrature for Sommerfeld contour:
    %        't' tanh-sinh curve, 'u' uniform, 'a' alpert real k [to implement] 
      if nargin<3, opts = []; end
      if isfield(opts,'omega'), ba.om = opts.omega;
      elseif isempty(ba.om), ba.om = ba.k; end               % default
      if isfield(opts,'nearsing'), ba.si = opts.nearsing;
      elseif isempty(ba.si), ba.si = ba.om/2; end            % "
      if isfield(opts,'shift'), ba.shift = opts.shift;
      elseif isempty(ba.shift), ba.shift = 0; end
      if isfield(opts,'maxy'), ba.maxy = opts.maxy;
      elseif isempty(ba.maxy), ba.maxy = 1; end
      if isfield(opts,'minx'), ba.minx = opts.minx;
      elseif isempty(ba.minx), ba.minx = 1; end
      if isfield(opts,'quad'), ba.quad = opts.quad;
      elseif isempty(ba.quad), ba.quad = 'b'; end
      % estimate parameters of contour...
      L = sqrt((log(1e-16)/ba.minx)^2 + ba.om^2); % Re[k] st x-decay to e_mach
      if ba.quad=='b'  % Alex's quadr tanh curve, with sinh bunching............
        d = min(4,6/ba.maxy);       % tanh imaginary dist (scales as O(1/maxy))
        de = d;          % max(d,ba.si), width-scale of tanh: keep slope <=1
        % params of remapping through sinh...
        b = max(2,real(ba.si)^1.2); % b expansion growth rate
        g = min(real(ba.si)/2,1); % g bunching fac, for be->a var change params
        if isempty(N) || N==0,    % estimate N... om>>1 effect, om<<1 effect
          Ns = [100, 3*max(ba.om,3)*ba.maxy, 16*b*asinh(L/g/b)];
          N = ceil(max(Ns)); end
        h = b*asinh(L/g/b) / (N/2-1);  % stopping at beta s.t. alpha stops at L
        be = (-N/2:(N+1)/2)*h;         % beta values (periodic trap rule)
        be = be(1:N);                  % clip to be exactly N
        wj = [1/2 ones(1,numel(be)-2) 1/2]*h; % weights
        a = g*b*sinh(be/b); wj = wj.*g.*cosh(be/b);  % transform be->a
        % change variable from real axis to contour shape (w/ imag part)...
        ba.kj = a - 1i*d*tanh(a/de); ba.wj = wj.*(1-1i*d/de*sech(a/de).^2);
      elseif ba.quad=='u'  % uniform, crude ....................................
        if nargin<2 || isempty(N), N = 100; end                % default N
        h = L/(N/2-1); ba.kj = ((1:N)-N/2)*h;  % compound trapezoid rule
        ba.wj = ones(1,numel(ba.kj))*h;
      end
      if N<1, error('N (either given or computed) must be at least 1!'); end
      ba.N = N;
      ba.kj = ba.kj + ba.shift;    % cool contour x-shift, non-Sommerfeld
    end
    
    function updateN(b,N) % ................ overloads from basis
    % UPDATEN - Change basis set # quadr pts in proportion to an overall N.
      b.N = ceil(N * b.nmultiplier);
      if isempty(b.doms) || isnan(b.doms(1).k)
        error('ftylayerpot must affect a domain with known wavenumber!'); end
      b.requadrature(b.N, struct('omega', b.doms(1).k)); % use domain's omega
    end
    
    function [A Ax Ay] = eval(b, p, o) % .........basis evaluator at points p
% EVAL - evaluate layer potential on pointset or on a segment, with jump rel
%
%  A = EVAL(bas, p, opts) returns a matrix mapping degrees of freedom in the
%   discretization of the Sommerfeld integral to the values on p the target
%   pointset or segment.
%
%  [A An] = EVAL(bas, p, opts) also returns normal derivatives using normals
%   in p.
%
%  [A Ax Ay] = EVAL(bas, p, opts) instead returns x- and y-derivatives,
%   ignoring the normals in p.
%
%  Options: opts.side = col vec of signs overriding which side eval occurs on
%               (+1, -1 means eval pt is to right, resp. left, of FTyLP).
%
% TODO: * make it use J-filter, or a version which returns J coeffs instead.
%         (would need to specify xloc of origin for J exp).
      if nargin<3, o = []; end
      om = b.k;               % method gets frequency om from affected domain
      N = numel(b.kj);                                   % # quadr pts
      M = numel(p.x);                                    % # target pts
      d = p.x(:) - b.orig;                               % col vec of targ rels
      x = real(d); y = imag(d);                          % eval pts rel to orig
      som = sqrt(om^2 - b.kj.^2);                        % Sommerfeld, row vec
      Aexp = exp(1i*(y*b.kj + abs(x)*som));  % matrix of exp factor in integrand
      if isfield(o, 'side'), signx = o.side; else signx = sign(x); end
      signx(find(signx==0)) = +1;                        % break symm if needed
      xnonneg = isempty(find(signx==-1)); xnonpos = isempty(find(signx==+1));
      
      if xnonneg
        A =  Aexp .* repmat((1i*b.a(1)./som + b.a(2)).*b.wj/2, size(x));
      elseif xnonpos
        A =  Aexp .* repmat((1i*b.a(1)./som - b.a(2)).*b.wj/2, size(x));
      else                                   % general mixed signs of x case
        A = Aexp .* (repmat(1i*b.a(1)./som.*b.wj/2,size(x)) + ...
                     (b.a(2)*signx/2)*b.wj);     % make faster when b.a(2)=0?
      end
      if nargout>1        % any derivs needed? if so, compute x-,y-derivs
        if xnonneg
          Ax = Aexp .* repmat((-b.a(1) + 1i*b.a(2).*som).*b.wj/2, size(x));
        elseif xnonpos
          Ax = Aexp .* repmat((b.a(1) + 1i*b.a(2).*som).*b.wj/2, size(x));
        else
          Ax = Aexp .* (-(b.a(1)*signx/2)*b.wj + ...
                        repmat(1i*b.a(2).*som.*b.wj/2, size(x)));
        end
        Ay = A .* repmat(1i*b.kj, size(x));   % y-deriv simple mult in k_y
      end
      if nargout==2                   % get An normal derivs from x-,y-derivs
        Ax = repmat(real(p.nx(:)), size(b.wj)).*Ax + ...
             repmat(imag(p.nx(:)), size(b.wj)).*Ay;
      end
    end % func
    
    function Q = evalftystripdiscrep(b, st) % ........... overloads from basis
    % EVALFTYSTRIPDISCREP - Q submatrix giving FTy LP effect on FTy discrepancy 
    %
    % Q = EVALFTYSTRIPDISCREP(bas, st). bas is the basis set, st the strip
    %   object. Currently supports standard ftylp's on L and copies on L+e, with
    %   buffer=0,1... Q is a 2NxN matrix, with NxN diagonal subblocks, since
    %   transl ops are diag in FTy rep.
    % 
    % ISSUES: check the non-x-parallel e case.
    %
    % * Obviously should make a multiply-by-Q routine since it's block diagonal!
      if ~isa(st, 'qpstrip'), error('st must be a qpstrip object!'); end
      N = b.Nf;
      om = b.k;                               % get omega from basis (ie domain)
      a = st.a;                               % Bloch phase
      d = real(st.e) * (1+2*st.buffer);       % x-displacement across strip
      som = sqrt(om^2 - b.kj.^2);             % Sommerfeld, row vec
      exf = exp(1i*som*d);                    % decay factors, row vec
      Q = [diag((b.a(1)*1i*(a-1/a)./som + b.a(2)*(-a-1/a)).*exf/2); ...
           diag((b.a(1)*(a+1/a) + b.a(2)*1i*(a-1/a)*som).*exf/2)];
      i = diagind(Q); Q(i) = Q(i) + b.a(2);   % jump relations add identities
      i = diagind(Q,N); Q(i) = Q(i) - b.a(1);
    end
    
    function showgeom(b, opts) % .................. crude show location of ray
      vline(real(b.orig), 'k-');
    end
    
  end % methods
  
  methods (Static) %......................
  end % methods
end