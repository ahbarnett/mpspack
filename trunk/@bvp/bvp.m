% BVP - define a Helmholtz or Laplace boundary value problem from domains
%
%  pr = bvp(doms) creates a BVP object pr which comprises the domains doms
%   (doms is a row vec of domain handles), using their segments and the
%   boundary conditions they carry, and the basis sets associated with the
%   domains.

classdef bvp < problem & handle
  properties
    rhs                             % right hand side (sqrt(w).f) col vector
  end
  
  methods % ---------------------------------- methods particular to BVPs
    function pr = bvp(doms) % .................... constructor
      if nargin==0, return; end    % needs empty constructor
      pr.doms = doms;
      pr.segs = [];             % now build a list of segments, crude O(N^2)
      for d=doms
        for s=d.seg
          if isempty(find(pr.segs==s))  % only add the segment if it's new
            pr.segs = [pr.segs s];
          end
        end
      end
    end

    function rhs = fillrighthandside(pr)   % ............... make RHS vector
      if isempty(pr.sqrtwei), error('must fill quadrature weights first'); end
      rhs = [];                 % get ready to stack stuff onto it as a big col
      for s=pr.segs
        if s.bcside==0          % matching condition
          if isnumeric(s.f)
            rhs = [rhs; s.f; s.g];           % data vec, stack as one big col
          else
            rhs = [rhs; s.f(s.t); s.g(s.t)]; % func of t, stack as one big col
          end
        elseif s.bcside==1 | s.bcside==-1    % BC (segment dofs natural order)
          if isnumeric(s.f)
            rhs = [rhs; s.f];       % data vector, stack as one big column
          else
            rhs = [rhs; s.f(s.t)];  % function of t, stack as one big column
          end
        end
      end
      rhs  = rhs .* pr.sqrtwei.';          % include sqrt quad weights factor
      if nargout==0, pr.rhs = rhs; end     % this only stores pr.rhs if no outp
    end
    
    function co = linsolve(pr) % ........................ linear solve
      co = pr.A \ pr.rhs;
      if nargout==0, pr.co = co; end   % only store internally if no output
    end
    
    function co = solvecoeffs(pr) % ................. fill and linear solve
      pr.fillquadwei;
      pr.fillrighthandside;
      pr.fillbcmatrix;
      co = pr.linsolve;
      if nargout==0, pr.co = co; end   % only store internally if no output
    end
    
    function r = bcresidualnorm(pr)    % ............... L2 bdry error
      r = norm(pr.A * pr.co - pr.rhs);
    end
    
  end
end
