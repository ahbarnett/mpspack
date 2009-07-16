classdef mfsbasis < handle & basis
% MFSBASIS - create a fundamental solutions (MFS) basis set
%
%  b = MFSBASIS(p, [], [], opts) where p is a pointset charge points
%   allows the user to choose the charge points, from which N is determined.
%
%  b = MFSBASIS(y, [], [], opts) where y is a row vector of charge points
%   allows the user to choose the charge points, from which N is determined.
%
%  b = MFSBASIS(Z, tau, N, opts) creates an MFS basis using charge points
%   at Z(t + i.tau) where t are equally distributed in [0,2pi], Z is a given
%   2pi-periodic analytic function, and tau is the distance off the curve.
%   N is the number of charge points, and k is wavenumber. opts may contain:
%   opts.real: if true use Y_0, otherwise use H_0, as fundamental solution.
%
%  b = MFSBASIS() creates MFS basis without any charge pts.
%
% Issues/notes:
%  * There should be other input formats for MFSBASIS!
%  * other location methods, inferface to mfsloc.m code from analytic_mfs
%  * Fully implement double and combined potential
%  * document opts.fast
%
% Also see: LAYERPOT
  properties
    Z                     % function handle for charge curve
    tau                   % imaginary distance
    t                     % row vec of parameter values in [0,2pi]
    y                     % row vec of charge points (C-#s)
    realflag              % if true, use real part only, ie Y_0
    fast                  % Hankel evaluation method (0=matlab; 1,2=faster)
    eta                   % If eta=inf (default) use single layer potential 
                          % MFS basis. For eta\neq inf
                          % use linear combination of double and single
                          % layer potential MFS basis
    %locmeth              % charge point location method
  end
  
  methods
    function b = mfsbasis(Z, tau, N, opts) % ............... constructor
      if nargin<4, opts = []; end
      if ~isfield(opts, 'real'), opts.real = 0; end
      b.realflag = opts.real;
      if ~isfield(opts, 'fast'), opts.fast = 2; end
      if ~isfield(opts, 'eta'), opts.eta = inf; end
      if opts.fast==2 && exist('@utils/greengardrokhlinhank106.mexglx')~=3
        opts.fast = 1;               % downgrade the speed if 106 not available
      end
      b.fast = opts.fast;
      b.eta  = opts.eta;
      if (b.fast>0 && b.eta<inf), 
          error('MFSERROR','Fast evaluation of combined basis not implemented');
      end
      if ~isempty(k), b.k = k; end
      if isnumeric(Z)                         % Z contains y list of MFS pts
        N = numel(Z);
        b.y = reshape(Z, [1 N]);
        b.N = N;
      elseif ~isempty(Z)                      % Z is a function handle
        b.N = N;
        b.Z = Z; b.tau = tau;
        b.t = 2*pi*(1:N)/N;                   % create row vec
        b.y = Z(b.t + 1i*tau);
      end
    end
    
    function Nf = Nf(b)
      Nf = b.N;
    end
    
    function [A B C] = eval(b, p, opts)
    % EVAL - evaluate MFS basis values, and maybe derivatives, at set of points
    %
    % Code adapted from /home/alex/bdry/inclus/evalbasis.m
      derivs = nargout-1;              % derivs = 1 direc, derivs = 2 both x & y
      N = b.N; M = numel(p.x); k = b.k;
      d = repmat(p.x, [1 N]) - repmat(b.y, [M 1]); % displacement matrix, C#s
      r = abs(d);
      if b.fast   % MEX code; bypass the split into fundsol & fundsol_deriv
        if b.fast==1
          [A radderivs] = utils.greengardrokhlinhank103(k*r); % get H0(kr), H1
        elseif b.fast==2
          [A radderivs] = utils.greengardrokhlinhank106(k*r); % get H0, H1
        end
        A = (1i/4) * A;
        radderivs = (1i*k/4) * radderivs;
        % now copy code from below b.fast=0, except reusing radderivs...
        if derivs==1                     % want directional deriv
          ny = repmat(p.nx, [1 N]);      % identical cols given by bdry normals
          cosphi = real(conj(ny).*d) ./ r;    % dot prod <normal, displacement>
          clear d ny
          B = utils.fundsol_deriv(r, -cosphi, k, radderivs); % (target)-normal deriv, NB sign
          clear cosphi
          if b.realflag, B = real(B); end         % keeps only the Y-0 part
        elseif derivs==2                          % want x,y-derivs
          dor = d./r;                             % contains x as re, y as im
          clear d
          [B radderivs] = utils.fundsol_deriv(r, -real(dor), k, radderivs); % x deriv 
          C = utils.fundsol_deriv(r, -imag(dor), k, radderivs); % y reuse rad part
          clear radderivs dor
          if b.realflag, B = real(B); C = real(C); end % keeps only the Y-0 part
        end
        
      else         % original, use matlab hankels in fundsol & fundsol_deriv...
      if derivs==1                     % want directional deriv
        ny = repmat(p.nx, [1 N]);      % identical cols given by bdry normals
        cosphi = real(conj(ny).*d) ./ r;      % dot prod <normal, displacement>
        clear d ny
        B = utils.fundsol_deriv(r, -cosphi, k); % (target)-normal deriv, NB sign
        clear cosphi
        if b.realflag, B = real(B); end         % keeps only the Y-0 part
      elseif derivs==2                          % want x,y-derivs
        dor = d./r;                             % contains x as re, y as im
        clear d
        [B radderivs] = utils.fundsol_deriv(r, -real(dor), k); % x deriv 
        C = utils.fundsol_deriv(r, -imag(dor), k, radderivs); % y reuse rad part
        clear radderivs dor
        if b.realflag, B = real(B); C = real(C); end % keeps only the Y-0 part
      end
      A = utils.fundsol(r, k);
      end
      if b.realflag, A = real(A); end           % keeps only the Y-0 part    
    end
    
    function showgeom(bas, opts) % .................. crude show MFS pts, etc
      if nargin<2, opts = []; end
      plot(real(bas.y), imag(bas.y), 'r+');
      if isfield(opts, 'label')
        n = ceil(numel(bas.y)/2);
        text(real(bas.y(n)), imag(bas.y(n)), opts.label);
      end
    end % func
        
  end % methods
end
