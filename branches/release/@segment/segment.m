% SEGMENT - create segment object
%
%  s = SEGMENT(M, [xi xf]) creates a line segment object from xi to xf, both
%  complex numbers.
%
%  s = SEGMENT(M, [xc R ti tf]) creates a circular arc segment with center
%  xi, radius R and angles from ti to tf. The order is important. If tf>ti
%  the orientation is counter-clockwise, otherwise clockwise.
%  
%  s = SEGMENT(M, {Z, Zp}) creats an analytic curve given by the image of 
%  the analytic function Z:[0,1]->C. Zp must be the derivative of Z. 
%  Z and Zp are function handles. 
%  
%  s = SEGMENT(M, {Z, Zp, Zpp}) works as above but also takes the second
%  derivative Z'' of Z. This is useful for layer potentials.
% 
%  s = SEGMENT(M, p, qtype) where p is any of the above, chooses quadrature type
%   qtype = 'p': periodic trapezoid (appropriate for periodic segments, M pts)
%           't': trapezoid rule (ie, half each endpoint, M+1 pts)
%           'c': Clenshaw-Curtis (includes endpoints, M+1 pts)
%           'g': Gauss (takes O(M^3) to compute, M pts)
%
%  If M is empty, a default value of 20 is used.
%
%  s = SEGMENT() creates an empty segment object.
%
% See also: POINTSET, segment/PLOT
%
% Copyright (C) 2008, 2009, Timo Betcke, Alex Barnett


classdef segment < handle & pointset
    properties
        t                      % quadrature parameter values in [0,1] (col vec)
        w                      % quadrature weights, sum = seg len (row vec)
        speed                  % |dZ/dt| at quadrature pts (col vec)
        kappa                  % kappa curvature at quad pts (col vec; optional)
        relaccel               % tangential accel, Maue-Kress (col vec; opt)
        qtype                  % quadrature type ('c', 'p', etc)
        eloc                   % [start point; end point] as C-#s
        eang                   % [start angle; end angle] as C-#s on unit circle
        Z                      % analytic function handle Z(t) on [0,1]
        Zp                     % derivative function handle dZ/dt on [0,1]
        Zpp                    % 2nd deriv d^2Z/dt^2 on [0,1] (optional)
        Zn                     % unit normal function handle on [0,1]
        approxv                % vertex list for polygonal approximation
        dom                    % domain handles bordered on + & - sides (2-cell)
        bcside                 % side +1/-1 BC is on, 0 for matching, or NaN
        a                      % BC value coeff (1-by-1), or +/- sides (1-by-2)
        b                      % BC n-deriv coeff (1-by-1), or +/- (1-by-2)
        f, g                   % BC data funcs or samples (f=value, g=n-deriv)
   end
    methods
      function s = segment(M, p, qtype)
        if nargin==0, return; end               % empty constructor (for copy)
        if nargin<3, qtype='c'; end             % default quadrature type
      
        % convert different types of input format all to an analytic curve...
        if iscell(p)         % ------------ analytic function (cell array)
          s.Z = p{1};        % use passed-in analytic func handles Z, Zp {,Zpp}
          s.Zp = p{2};
          if numel(p)>2
            s.Zpp = p{3};
          end
          Napprox = 100;     % # pts for crude inside-polygon test, must be even
        elseif numel(p)==2   % ------------ straight line
          d = p(2)-p(1);
          s.Z = @(t) p(1) + d*t;
          s.Zp = @(t) d + 0*t;     % constant (0*t trick to make size of t)
          s.Zpp = @(t) 0*t;
          Napprox = 1;
        elseif numel(p)==4   % ------------- arc of circle
          s.Z = @(t) p(1) + p(2)*exp(1i*(p(3) + t*(p(4)-p(3))));
          s.Zp = @(t) 1i*(p(4)-p(3))*p(2)*exp(1i*(p(3) + t*(p(4)-p(3))));
          s.Zpp = @(t) -(p(4)-p(3))^2*p(2)*exp(1i*(p(3) + t*(p(4)-p(3))));
          if abs(p(4)-p(3)-2*pi)<1e-15, qtype = 'p'; end % closed -> periodic q
          Napprox = 50;      % # pts for crude inside-polygon test, must be even
        else
          error('segment second argument not valid!');
        end
        s.qtype = qtype;             % keep for later reference
        s.requadrature(M, qtype);    % set up quadrature pts, w, nx, etc...
        s.Zn = @(t) -1i*s.Zp(t)./abs(s.Zp(t)); % new; supercedes normal method
        s.eloc = s.Z([0;1]);
        eZp = s.Zp([0;1]);                     % derivs at the 2 ends
        s.eang = eZp./abs(eZp);
        s.approxv = s.Z((0:Napprox)'/Napprox); % start, endpt (drop one later)
        s.dom = {[] []};             % not bordering any domains
        s.bcside = NaN;              % no BCs or matching conditions
      end
      
      function requadrature(segs, M, qtype)
      % REQUADRATURE - change a segment's quadrature scheme or number of points
      %
      %  REQUADRATURE(seg, M) changes the number of points to M (or M+1
      %   depending on the type, see below), without changing the quadrature
      %   type. If seg is a list of segments, it does so for each.
      %
      %  REQUADRATURE(seg, M, qtype) chooses new M and new quadrature type
      %   qtype = 'p': periodic trapezoid (appropriate for periodic segments,
      %                M pts are created)
      %           't': trapezoid rule (ie, half each endpoint, M+1 pts)
      %           'c': Clenshaw-Curtis (includes endpoints, M+1 pts)
      %           'g': Gauss (takes O(M^3) to compute, M pts)
      %
      %  If M is empty, a default value of 20 is used.
      %
      % See also: SEGMENT
        if isempty(M), M = 20; end              % default M
        for s=segs
          if nargin<3, qtype=s.qtype; end         % preserve quadrature type
          switch qtype % choose a quadrature rule function on [-1,1]...
           case 'p',
            quadrule = @quadr.peritrap;
           case 't',
            quadrule = @quadr.traprule;
           case 'c',
            quadrule = @quadr.clencurt;
           case 'g',
            quadrule = @quadr.gauss;  % note via eig returns increasing x order
           otherwise,
            error(sprintf('requadrature: unknown quadrature type %s!', qtype));
          end
          [z s.w] = quadrule(M); % NB must give monotonic increasing x in [-1,1]
          s.t = (1+z)/2;               % t in [0,1], increasing
          s.x = s.Z(s.t);
          dZdt = s.Zp(s.t);            % Z' eval at t
          s.speed = abs(dZdt);
          s.w = s.w/2 .* s.speed.';  % quadr weights wrt arclength on segment
          s.nx = -1i*dZdt./s.speed;
          if ~isempty(s.Zpp)         % d^2Z/dt^2 available...
            s.kappa = -real(conj(-1i*dZdt).*s.Zpp(s.t)) ./ s.speed.^3;%curvature
            s.relaccel = real(conj(dZdt).*s.Zpp(s.t)) ./ s.speed.^2; %Maue-Kress
          end
        end
      end
      
      function setbc(seg, pm, a, b, f) % .......................... setBC
      % SETBC - set a (in)homogeneous boundary condition on a segment
      %
      %  setbc(seg, pm, type) creates a homogeneous BC on the side given by pm
      %   (+1 is the postive normal side, -1 is the `interior domain' or back
      %   side; the sign of pm is opposite for domain constructor), of type 'd'
      %   or 'n' for Dirichlet or Neumann. Checks if a domain is attached on the
      %   requested side, reports error if not. seg and pm may also be arrays of
      %   segment handles and signs (if pm is length 1 it will be applied to
      %   each segment).
      %
      %  setbc(seg, pm, a, b) imposes the BC au + bu_n = 0 on the side pm, using
      %   normal sense as above. a, b are constants (may be functions in future)
      %
      %  setbc(seg, pm, type, [], f) imposes inhomogeneous BC of type 'd' or 'n'
      %   with data func f(t) where t in [0,1] parametrizes the segment, if f is
      %   a function handle with one argument. If f is a function handle with
      %   two arguments, it is interpreted as f(x,y) for points (x,y) on the
      %   segment. If instead f is a (col vec) array of same size as
      %   seg.x, is interpreted as data sampled on the quadrature points.
      % 
      %  Notes/Issues:
      %  1) a seg can currently carry only one BC, or instead a matching
      %   condition (which is a pair of relations on the segment).
      %  2) this routine merely copies info into the segment. The chief routines
      %   which interpret and calculate using the BC/matching include
      %   problem.fillbcmatrix, bvp.fillrighthandside, and
      %   scattering.setincidentwave.
      %  3) Having f(x,y) vs f(t) distinguish whether f is a func of location
      %   or of parameter t is not great. Prefer f(z) vz f(t) - how distinguish?
        if numel(pm)==1, pm = pm*ones(size(seg)); end
        for i=1:numel(seg)
          s = seg(i);
          ind = (1-pm(i))/2+1;            % index 1 (+) or 2 (-)
          if isempty(s.dom(ind))
            error(sprintf('side %d of seg not connected to a domain!', pm(i)));
          end
          s.bcside = pm(i);               % same natural segment side convention
          if nargin<5, s.f = @(t) zeros(size(s.t)); % covers homogeneous case
          elseif isa(f,'function_handle')
            s.f = f;                      % func handle (not eval'd until need)
            try, f(0); catch me,          % test to see if f has >1 argument
              if strcmp(me.identifier, 'MATLAB:inputArgUndefined') % hack
                s.f = @(t) f(real(s.Z(t)),imag(s.Z(t))); % interpret as f(x,y)
              end
            end
          else, s.f = f;                  % will be an array
          end
          if nargin<4 | isempty(b)
            switch a
             case {'D', 'd'}
              s.a = 1; s.b = 0;
             case {'N', 'n'}
              s.a = 0; s.b = 1;
             otherwise
              error(sprintf('unknown BC type: %s', a))
            end
          else
            s.a = a; s.b = b;
          end
        end
      end % func
      
      function setmatch(seg, a, b, f, g) % ....................... setmatch
      % SETMATCH - set (in)homogeneous matching conditions across a segment
      %
      %  setmatch(seg, a, b) creates a homogeneous matching condition
      %   relating values and normal derivatives on the two sides of a segment.
      %   a, b are 1-by-2 arrays [a^+ a^-] and [b^+ b^-].
      %   (+ is natural positive normal of segment; - is 'interior' back side).
      %   If seg is an array of segments the same matching is applied to each.
      %   The two conditions imposed are
      %                                   a^+ u^+ + a^- u^-     = 0
      %                                   b^+ u_n^+ + b^- u_n^- = 0
      %
      %  setmatch(seg, a, b, f, g) replaces the right hand sides of the
      %   above by inhomogeneous matching functions f,g of t in [0,1]
      %
      %  setmatch(seg, 'diel', pol) uses the dielectric constants of the
      %   bordering domains to set (a,b) as above appropriate for a dielectric
      %   matching condition. pol is either 'tm' or 'te' for
      %   transverse-magnetic (E_z is scalar) or transverse-electric (H_z is
      %   scalar). f and g may be appended as above. 
      %
      % Notes: see Notes for SETBC
      %
      % See also DIELECTRICCOEFFS.
        for i=1:numel(seg)
          s = seg(i);
          if isempty(s.dom{1}) | isempty(s.dom{2})
            error('both sides of seg must be connected to a domain!');
          end
          s.bcside = 0;
          if nargin<4, f = @(t) zeros(size(s.t)); end    % covers homog f case
          if nargin<5, g = @(t) zeros(size(s.t)); end    % covers homog g case
          s.f = f; s.g = g;
          if strcmp(a,'diel')                 % use dielectrics to set matching
            [s.a s.b] = segment.dielectriccoeffs(b, s.dom{1}.refr_ind, ...
                                                 s.dom{2}.refr_ind);
          else                                % use a,b passed in, override diel
            if numel(a)~=2 | numel(b)~=2, error('a and b must be 1-by-2!'); end
            s.a = a; s.b = b;
          end
        end
      end % func
      
      function newseg = scale(seg, fac)  % ................... dilate a segment
      % SCALE - rescale (dilate) a segment (or list) about the origin
      %
      %  scale(seg, fac) changes the segment seg (or an array of segments)
      %   to be rescaled about the origin by factor fac>0.
      %
      %  newseg = scale(seg, fac) makes a new segment (or list) which is seg
      %   rescaled about the origin by factor fac>0.
        if nargout>0                             % make a duplicate
          newseg = [];
          for s=seg; newseg = [newseg utils.copy(s)]; end
        else
          newseg = seg;                          % copy handle, modify original
        end
        for s=newseg
          s.x = fac * s.x;
          s.w = fac * s.w;
          s.speed = fac * s.speed;
          s.kappa = s.kappa / fac;
          Z = s.Z; Zp = s.Zp;
          s.Z = @(t) fac * Z(t);
          s.Zp = @(t) fac * Zp(t);
          if ~isempty(s.Zpp)
            Zpp = s.Zpp;
            s.Zpp = @(t) fac * Zpp(t);
          end
           s.eloc = fac * s.eloc;
          s.approxv = fac * s.approxv;
        end
      end % func
      
      function newseg = translate(seg, a)  % ............. translate a segment
      % TRANSLATE - translate a segment (or list of segments)
      %
      %  translate(seg, a) changes the segment (or array of segments) seg by
      %   translating by the complex number a.
      %
      %  newseg = translate(seg, a) instead returns a new segment (or list)
      %   obtained by translating the segment (or list) seg.
      %
      %  Note that multiple translations cause recursive depth in seg.Z, slow!
        if nargout>0                             % make a duplicate
          newseg = [];
          for s=seg; newseg = [newseg utils.copy(s)]; end
        else
          newseg = seg;                          % copy handle, modify original
        end
        for s=newseg
          s.x = s.x + a;
          Z = s.Z;
          s.Z = @(t) Z(t) + a;
          s.eloc = s.eloc + a;
          s.approxv = s.approxv + a;
        end
      end % func
      
     function newseg = rotate(seg, t)  % ............. rotate a segment
      % ROTATE - rotate a segment (or list of segments) about the origin
      %
      %  rotate(seg, t) changes the segment (or array of segments) seg by
      %   rotating CCW by the angle t
      %
      %  newseg = rotate(seg, t) instead returns a new segment (or list)
      %   obtained by rotating the segment (or list) seg.
        a = exp(1i*t);                           % a is on unit circle
        if nargout>0                             % make a duplicate
          newseg = [];
          for s=seg; newseg = [newseg utils.copy(s)]; end
        else
          newseg = seg;                          % copy handle, modify original
        end
        for s=newseg
          s.x = a * s.x; s.nx = a * s.nx;
          Z = s.Z; Zp = s.Zp; Zn = s.Zn;
          s.Z = @(t) a * Z(t); s.Zp = @(t) a * Zp(t); s.Zn = @(t) a * Zn(t);
          if ~isempty(s.Zpp)
            Zpp = s.Zpp;
            s.Zpp = @(t) a * Zpp(t);
          end
          s.eloc = a * s.eloc; s.eang = a * s.eang;
          s.approxv = a * s.approxv;
        end
      end % func
      
     function newseg = reflect(seg, ax)  % ............. reflect a segment
      % REFLECT - reflect a segment (or list of segments) about x or y axis
      %
      % reflect(seg, ax) changes the segment (or array of segments) seg by
      %   reflecting through either ax='x' (default if ax empty or not given)
      %   or 'y' axis.
      %
      % newseg = reflect(seg, ax) instead returns a new segment (or list)
      %   obtained by reflecting the segment (or list) seg.
        if nargin<2 || isempty(ax), ax='x'; end              % default
        if nargout>0                             % make a duplicate
          newseg = [];
          for s=seg; newseg = [newseg utils.copy(s)]; end
        else
          newseg = seg;                          % copy handle, modify original
        end
        for s=newseg
          if ax=='x'
          s.x = conj(s.x); s.nx = -conj(s.nx);    % note sense change!
          Z = s.Z; Zp = s.Zp; Zn = s.Zn;
          s.Z = @(t) conj(Z(t));
          s.Zp = @(t) conj(Zp(t)); s.Zn = @(t) -conj(Zn(t));  % note sense!
          if ~isempty(s.Zpp)
            Zpp = s.Zpp;
            s.Zpp = @(t) conj(Zpp(t));
          end
          s.eloc = conj(s.eloc); s.eang = conj(s.eang);
          s.approxv = conj(s.approxv);
          elseif ax=='y'
            error('y-reflection not yet implemented');
          else
            error('unknown reflection axis!');
          end
        end
      end % func
      
     function disconnect(segs)  % ...... disconnect a seg list from any domains
     % DISCONNECT - disconnect a segment or segment list from any domains
        for s=segs
          s.dom = {[] []};             % not bordering any domains
        end
      end

      h = plot(s, pm, o)
    end % methods

    % --------------------------------------------------------------------
    methods(Static)    % these don't need segment obj to exist to call them...
      s = polyseglist(M, p, qtype)
      s = radialfunc(M, fs)
      s = smoothstar(M, a, w)
      [a b] = dielectriccoeffs(pol, np, nm)
    end % methods
end

