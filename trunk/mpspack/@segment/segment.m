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
        w                      % quadrature weights, sum = seg len (row vec)
        eloc                   % [start point; end point] as C-#s
        eang                   % [start angle; end angle] as C-#s on unit circle
        Z                      % analytic function Z(t) on [0,1]
        Zp                     % derivative dZ/dt on [0,1]
        approxv                % vertex list for polygonal approximation
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
        t = (1+z)/2;            % t in [0,1], increasing
        s.x = s.Z(t);
        dZdt = s.Zp(t);         % Z' eval at t
        s.w = s.w/2 .* abs(dZdt).';  % quadr weights wrt arclength on segment
        s.nx = -1i*dZdt./abs(dZdt);
        s.eloc = s.Z([0;1]);
        eZp = s.Zp([0;1]);                     % derivs at the 2 ends
        s.eang = eZp./abs(eZp);
        s.approxv = s.Z((0:Napprox)'/Napprox); % start, endpt (drop one later)
      end
      
      function n = normal(s, t) % ................. returns normal at t in [0,1]
        dZdt = s.Zp(t);
        n = -1i*dZdt./abs(dZdt);
      end
      
      h = plot(s, pm, o)
    end

    % --------------------------------------------------------------------
    methods(Static)    % these don't need segment obj to exist to call them...
      s = polyseglist(M, p)
    end
end

