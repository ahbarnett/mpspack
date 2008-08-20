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
  
    function [Q data] = evalunitcelldiscrep(b, uc, opts)
    % EVALUNITCELLDISCREP - return Q matrix which maps basis coeffs to discrep
    %
    %  Default routine for bases which don't know anything special about the UC
    %
    %  Q = EVALUNITCELLDISCREP(bas, uc, opts) uses options in opts including:
    %   opts.dom: the domain in which the basis is evaluated (eg, conn. dom)
    %   opts.wei = 0,1,2...: 0 has no source nei copies, 1 has 3x3, etc.
    %
    % Notes: 1) now is simplified, since no layerpotside needed.
    %        2) To do: reuse if opts.data passed in, or make uc.Cdata storage!
      if ~isa(uc, 'qpunitcell'), error('uc must be a qpunitcell object!'); end
      savemats = (nargout>1);
      if nargin<3, opts = []; end
      if ~isfield(opts, 'nei'), opts.nei = 0; end
      nei = opts.nei;
      opts.dom = uc;
      if nei==0                           % no neighboring local copies
        [L Ln] = b.eval(uc.L, opts);
        [R Rn] = b.eval(uc.R, opts);
        Q = [L - (1/uc.a)*R; Ln - (1/uc.a)*Rn];           % formula for discrep
        if savemats, data.L = L; data.Ln = Ln; data.R = R; data.Rn = Rn; end
        clear L Ln R Rn
        [B Bn] = b.eval(uc.B, opts);
        [T Tn] = b.eval(uc.T, opts);
        Q = [Q; B - T*(1/uc.b); Bn - Tn*(1/uc.b)];       % formula for discrep
        if savemats, data.B = B; data.Bn = Bn; data.T = T; data.Tn = Tn; end
      else                                % sum local neighbors w/ phases
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
    
  end % methods
end
