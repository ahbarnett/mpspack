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
  
    function [Q data] = evalunitcelldiscrep(b, uc, opts)
    % EVALUNITCELLDISCREP - return Q matrix which maps basis coeffs to discrep
    %
    %  Default routine for bases which don't know anything special about UC
    %
    % * to do: reuse if opts.data!
      if ~isa(uc, 'qpunitcell'), error('uc must be a qpunitcell object!'); end
      savemats = (nargout>1);
      opts.layerpotside = 1; [L Ln] = b.eval(uc.L, opts); % LP signs: uc's -pm
      opts.layerpotside = -1; [R Rn] = b.eval(uc.R, opts);
      Q = [L - (1/uc.a)*R; Ln - (1/uc.a)*Rn];             % formula for discrep
      if savemats, data.L = L; data.Ln = Ln; data.R = R; data.Rn = Rn; end
      clear L Ln R Rn
      opts.layerpotside = 1; [B Bn] = b.eval(uc.B, opts);
      opts.layerpotside = -1; [T Tn] = b.eval(uc.T, opts);
      Q = [Q; B - T*(1/uc.b); Bn - Tn*(1/uc.b)];          % formula for discrep
      if savemats, data.B = B; data.Bn = Bn; data.T = T; data.Tn = Tn; end
    end % func
    
  end % methods
end
