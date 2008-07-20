% SEGMENT - create segment object
%
%  s = SEGMENT(M, [xi xf]) create line segment object from xi to xf, both C#s.
%
%  s = SEGMENT(M, [xc R ti tf]) create circular arc segment, center xc (C-#),
%   radius R, from angle ti to tf. Order is important: if tf>ti then goes CCW,
%   otherwise CW.
%
%  s = SEGMENT(M, {Z, Zp}) analytic curve given by image of analytic function
%   Z:[0,1]->C. Zp must be the derivative function Z'. Note the argument is a
%   1-by-2 cell array of function handles.
%
%  s = SEGMENT(M, p, qtype) where p is any of the above, chooses quadrature type
%   qtype = 'p': periodic trapezoid (appropriate for periodic segments, M pts)
%           't': trapezoid rule (ie, half each endpoint, M+1 pts)
%           'c': Clenshaw-Curtis (includes endpoints, M+1 pts)
%           'g': Gauss (takes O(M^3) to compute, M pts)
%
% See also: POINTSET, segment/PLOT

classdef segment < handle & pointset
    properties
        t                      % quadrature parameter values in [0,1]
        w                      % quadrature weights, sum = seg len (row vec)
        eloc                   % [start point; end point] as C-#s
        eang                   % [start angle; end angle] as C-#s on unit circle
        Z                      % analytic function handle Z(t) on [0,1]
        Zp                     % derivative function handle dZ/dt on [0,1]
        Zn                     % unit normal function handle on [0,1]
        approxv                % vertex list for polygonal approximation
        dom                    % domain handles bordered on + & - sides (2-cell)
        domseg                 % which segment # in bordering domains (1-by-2)
        bcside                 % side +1/-1 BC is on, 0 for matching, or NaN
        a                      % BC value coeff (1-by-1), or +/- sides (1-by-2)
        b                      % BC n-deriv coeff (1-by-1), or +/- (1-by-2)
        f, g                   % BC data funcs or samples (f=value, g=n-deriv)
   end
    methods
      function s = segment(M, p, qtype)
        if nargin<3, qtype='c'; end             % default quadrature type
      
        % convert different types of input format all to an analytic curve...
        if iscell(p)         % ------------ analytic function (cell array)
          s.Z = p{1};        % use passed-in analytic func handles Z, Zp
          s.Zp = p{2};
          Napprox = 100;     % # pts for crude inside-polygon test, must be even
        elseif numel(p)==2   % ------------ straight line
          d = p(2)-p(1);
          s.Z = @(t) p(1) + d*t;
          s.Zp = @(t) d + 0*t;     % constant (0*t trick to make size of t) 
          Napprox = 1;
        elseif numel(p)==4   % ------------- arc of circle
          s.Z = @(t) p(1) + p(2)*exp(1i*(p(3) + t*(p(4)-p(3))));
          s.Zp = @(t) 1i*(p(4)-p(3))*p(2)*exp(1i*(p(3) + t*(p(4)-p(3))));
          Napprox = 50;      % # pts for crude inside-polygon test, must be even
        else
          error('segment second argument not valid!');
        end
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
          error(sprintf('segment: unknown quadrature type %s!', qtype));
        end
        [z s.w] = quadrule(M);  % NB must give monotonic increasing x in [-1,1]
        s.t = (1+z)/2;            % t in [0,1], increasing
        s.x = s.Z(s.t);
        dZdt = s.Zp(s.t);         % Z' eval at t
        s.w = s.w/2 .* abs(dZdt).';  % quadr weights wrt arclength on segment
        s.nx = -1i*dZdt./abs(dZdt);
        s.Zn = @(t) -1i*s.Zp(t)./abs(s.Zp(t)); % new; supercedes normal method
        s.eloc = s.Z([0;1]);
        eZp = s.Zp([0;1]);                     % derivs at the 2 ends
        s.eang = eZp./abs(eZp);
        s.approxv = s.Z((0:Napprox)'/Napprox); % start, endpt (drop one later)
        s.dom = {[] []};             % not bordering any domains
        s.domseg = [0 0];            % dummy segment #s in bordering domains
        s.bcside = NaN;              % no BCs or matching conditions
      end
      
      function n = normal(s, t) % ................. returns normal at t in [0,1]
      % NORMAL - compute unit normal as C-# given segment parameter t in [0,1]
      %
      %  Note: this is duplicated by the property Zn in segment object
        dZdt = s.Zp(t);
        n = -1i*dZdt./abs(dZdt);
      end
      
      function setbc(s, pm, a, b, f) % .......................... setBC
      % SETBC - set a (in)homogeneous boundary condition for a segment
      %
      %  setbc(seg, pm, type) creates a homogeneous BC on the side given by pm
      %   (+1 is inside, ie opposite from normal; -1 is other side), of type 'd'
      %   or 'n' for Dirichlet or Neumann. Normal sense for pm=-1 is same as
      %   the segment natural sense. Checks if a domain is attached on the
      %   requested side, reports error if not.
      %
      %  setbc(seg, pm, a, b) imposes the BC au + bu_n = 0 on the side pm, using
      %   normal sense as above. a, b are constants (may be functions in future)
      %
      %  setbc(seg, pm, type, [], f) imposes inhomogeneous BC of type 'd' or 'n'
      %   with data func f(t) where t in [0,1] parametrizes the segment, if f is
      %   a function handle. If instead f is a (col vec) array of same size as
      %   seg.x, is interpreted as data sampled on the quadrature points.
      %
      %  Note that a seg can currently carry only one BC, or instead a matching
      %   condition (which is a pair of relations on the segment).
        ind = (1-pm)/2+1;            % index
        if isempty(s.dom(ind))
          error(sprintf('side %d of seg not connected to a domain!', pm));
        end
        s.bcside = pm;               % should be +-1
        if nargin<5, f = @(t) zeros(size(s.t)); end   % covers homogeneous case
        s.f = f;           % may be func handle (not eval'd until need) or array
        if isempty(b)
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
      end % func
      
      function setmatch(s, a, b, f, g) % ....................... setmatch
      % SETMATCH - set (in)homogeneous matching conditions across a segment
      %
      %  setmatch(seg, a, b) creates a homogeneous matching condition
      %   relating values and normal derivatives on the two sides of a segment.
      %   a, b are 1-by-2 arrays [a^+ a^-] and [b^+ b^-].
      %   (+ is inside, ie opposite from normal; - is other side)
      %   The two conditions imposed are
      %                                   a^+ u^+ + a^- u^-     = 0
      %                                   b^+ u_n^+ + b^- u_n^- = 0
      %
      %  setmatch(seg, a, b, f, g) replaces the right hand sides of the
      %   above by inhomogeneous matching functions f,g of t in [0,1]
        if isempty(s.dom(1)) | isempty(s.dom(2))
          error('both sides of seg must be connected to a domain!');
        end
        s.bcside = 0;
        if nargin<4, f = @(t) zeros(size(s.t)); end    % covers homog f case
        if nargin<5, g = @(t) zeros(size(s.t)); end    % covers homog g case
        s.f = f; s.g = g;
        if numel(a)~=2 | numel(b)~=2, error('a and b must be 1-by-2!');
          s.a = a; s.b = b;
        end
      end % func
      
      h = plot(s, pm, o)
    end % methods

    % --------------------------------------------------------------------
    methods(Static)    % these don't need segment obj to exist to call them...
      s = polyseglist(M, p)
      
      function disconnect(segs)  % ...... disconnect a seg list from any domains
        for s=segs
          s.dom = {[] []};             % not bordering any domains
          s.domseg = [0 0];            % dummy segment #s in bordering domains
        end
      end
    end % methods
end

