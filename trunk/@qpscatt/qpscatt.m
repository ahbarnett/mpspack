% QPSCATT - define a quasi-periodic frequency-domain scattering problem
%
%  pr = qpscatt(airdoms, doms, str) creates a scattering problem object pr
%   as with scattering, except a quasi-periodic version with periodicity of
%   the vertical strip domain str. The solution method will be via FTyLPs.
%
%  pr = qpscatt(airdoms, doms, str, opts) sets certain algorithm parameters:
%    opts.M = # degrees of freedom in each FTyLP
%    opts.nei = 0,1,2... how many direct sums to include on each side
%    opts.buf = 0,1,2... how wide the qpstrip domain is including buffer size
%
% Issues: * also have qpscatt create the qpstrip (better for buf>0) from args?
%
% See also PROBLEM, BVP, SCATTERING

classdef qpscatt < scattering & handle
  properties
    t                                       % the qpstrip object (stores alpha)
    nei
    buf
    M                                       % # nodes FTyLP Sommerfeld contour
    kpx, kpy                                % k-plane poles (k_x,k_y) (QP waves)
    kpyfix                                  % list of k_y for poles to cross
  end
  
  methods % -------------------------------- methods particular to qpscatt
    function pr = qpscatt(airdoms, doms, t, o) % ................... constructor
      pr = pr@scattering(airdoms, doms);
      if ~isempty([airdoms doms])              % if nonempty constructor
        pr.t = t;
        if nargin<4, o = []; end               % process method options
        if isfield(o,'M'), pr.M = o.M; else pr.M = 120; end  % defaults
        if isfield(o,'nei'), pr.nei = o.nei; else pr.nei = 2; end
        if isfield(o,'M'), pr.buf = o.buf; else pr.buf = 0; end
        % choose a periodizing basis for the qpstrip...
        pr.t.addqpftylayerpots(struct('M',pr.M,'nearsing',3)); % NB default!
      end
    end

    function setincidentwave(pr, t)
    % SETINCIDENTWAVE - set inc. plane wave direction, Bloch, k-fix, QP bases
    %
    %  setincidentwave(pr, t) sets up the Bloch alpha and BCs required in
    %   qpscatt problem pr, for planewave at angle t in [-pi,0].
    %   Wavenumber must already be set in the problem. FTyLPs are adjusted
    %   based on nearsing, and Wood's correction bases are chosen and set up.
    %
    %  Careful: calling this routine overwrites all inhomogeneity functions or
    %   data f, g stored on any of the problem's segments. However it preserves
    %   existing a, b BC or matching coeffs (which must be set up on entry).
      if nargin==2
        if isempty(pr.k), error('qpscatt problem needs wavenumber set'); end
        setincidentwave@scattering(pr, t); % call superclass method
        om = pr.k; kvec = om*exp(1i*t);
        pr.t.setbloch(exp(1i*real(conj(kvec) * pr.t.e))); % set Bloch alpha
        % set the k-poles list in pr object...
        d = abs(pr.t.e); n = 2*ceil(om*d/2/pi)+10;  % # diffraction orders?
        kpx = om*cos(t)+(-n:n)*2*pi/d; kpy = sqrt(om^2-kpx.^2); % kp is y-compt
        [dum i] = sort(abs(kpy),'ascend');       % note +ve sqrt key here
        pr.kpx = kpx(i); pr.kpy = kpy(i);
        % reset the QP FTyLP bases, and add any Wood's anomaly fix bases...
        for i=1:2, pr.t.bas{i}.requadrature(pr.t.bas{i}.N, ...
                      struct('omega',om, 'nearsing', abs(pr.kpy(1)))); end
        pr.t.bas = pr.t.bas(1:2); pr.kpyfix = []; % kill off any Wood fix bases
      end
    end
    
    function showkyplane(pr)
    % SHOWKYPLANE - plot Sommerfeld quadrature nodes, k_y poles, +-om, etc
       figure; plot(pr.t.bas{1}.kj, '+'); axis tight equal; hold on;
       plot([pr.kpy -pr.kpy], 'rx');
       plot(real(pr.kpyfix),imag(pr.kpyfix),'go');
       plot([-pr.k pr.k], [0 0], 'k*');
       title('complex k (ie k_y) plane: contour, poles'); 
    end
    
    function rhs = fillrighthandside(pr)
    % FILLRIGHTHANDSIDE - overloads BVP routine, setting mismatch and discrep
    %
    % Note: incident field is QP so leads to zero discrepancy.
      rhs = [];                 % get ready to stack stuff onto it as a big col
      obstrhs = fillrighthandside@bvp(pr); % obstacle mismatch part of RHS
      rhs = [obstrhs; zeros(pr.t.N,1)];    % QP block square, so use # QP dofs
      if nargout==0, pr.rhs = rhs; end     % this only stores pr.rhs if no outp
    end

    function A = fillbcmatrix(pr, opts)
    % FILLBCMATRIX - fills [A B;C Q] for QP scatt prob: dofs->mismatch/discrep
    %
    % Issues/Notes:
    %  * everything dense. Need replace Q by diagonal multiplication
    %    method, and make A a method not a matrix for FMM case.
    %    Better, make A a struct with mult method, B, C, and diag(Q) in it.
      bob = pr;                        % following fillbcmatrix@problem
      N = bob.setupbasisdofs;          % sets up obst dofs in pr
      t = pr.t;                        % the qpstrip domain
      Q = t.evalbasesdiscrep();        % is block-diag - need not really fill!
      % fill A: self, then directly sum neighbor images contribs...
      A = fillbcmatrix@problem(pr);    % obstdofs->mismatch block
      nei = pr.nei; buf = pr.buf; a = pr.t.a; % get Bloch alpha
      %for i=1:numel(bob.bas)                 % all bases in problem or basis obj
       % b = bob.bas{i}; ns = bob.basnoff(i)+(1:b.Nf); % dof indices for this bas
        
      %for n=[-nei:-1 1:nei], A = A + a^n * l.eval(pointset(s.x-d*n)); end

      if 0 % fill B: effect of FTyLPs on the airdom-touching segs in pr...
      for s=pr.segs
        if s.bcside==0                       % matching, may be dielectric
        if s.dom{1}.isair~=s.dom{2}.isair    % air-nonair junction
         % extract the relevant coeffs for whichever side have inc field on it:
          if s.dom{1}.isair       % + side has inc, - side has zero field
            a = s.a(1); b = s.b(1);
          else                    % - side has inc, + side has zero field
            a = s.a(2); b = s.b(2);
          end
          s.f = @(t) -a * ui(s.Z(t));        % jump in value of u_s (NB minus)
          s.g = @(t) -b * (uix(s.Z(t)).*real(s.Zn(t)) + ...
                uiy(s.Z(t)).*imag(s.Zn(t))); % jump in deriv of u_s
        end
        elseif s.bcside==1 | s.bcside==-1    % BC
          ind = (1-s.bcside)/2+1; % index 1 or 2 for which side the BC on
          d = s.dom{ind};         % handle of domain on the revelant side
          if d.isair              % air-to-metallic boundary
            if s.b==0             % Dirichlet only. Minus sign: u_s cancels u_i
              s.f = @(t) -s.a * ui(s.Z(t));  % NB keeping f a func not an array
            else                  % Robin, Neumann
              s.f = @(t) -s.a * ui(s.Z(t)) -s.b * (uix(s.Z(t)).*real(s.Zn(t))...
                              + uiy(s.Z(t)).*imag(s.Zn(t)));
            end
          else                    % internal-to-metallic boundary
            s.f = @(t) zeros(size(s.t));     % f = homogeneous BCs
          end
        end
      end % segs loop
end
      
      if 0
      
      % fill C....
      for i=1:numel(bob.bas)           % all bases in problem or basis obj
        b = bob.bas{i};
        ns = bob.basnoff(i)+(1:b.Nf);     % dof indices for this bas
% restrict to only those affecting airdoms
    
        % ...
      end
    
      if ~isempty(kpyfix), % augment matrix with Wood fix ...
        end
      
    end
    end
    
  end % methods
end
