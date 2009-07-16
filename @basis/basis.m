classdef basis < handle
    
  % Basis class - Abstract class that defines the interfaces which are
  % common for all basis objects.
  
  properties
    doms                    % array of handles of domains affected by basis
                            %     (added 7/14/09, required)
    N                       % Degree of basis set: interpretation depends on
                            %     type of basis (see basis.Nf method)
  end
  
  methods (Abstract)
    [A, A1, A2] = eval(b, pts) % Evaluate a basis on a set of points
                               % A=eval(pts) returns only fct. values
                               % [A A1]=eval(pts) returns fct. values & normal
                               % derivatives.
                               % [A A1 A2]=eval(pts) returns fct. values plus
                               % x and y derivatives.

    Nf = Nf(b, opts)           % method returning number of degrees of freedom
  end
  
  methods % ....................................... actual methods ...........
  
    function updateNf(b) % ................ dummy for all but layerpot bases
    end
    
    function k = k(b, opts) % .................... looks up the basis wavenumber
    % K - look up a basis set's wavenumber from the domain(s) it affects
    %
    % k = k(bas) returns the wavenumber in the domain affected by basis set bas.
    %
    % Issues: need to generalize to basis sets that affect >1 domain.
      k = b.doms(1).k;       % note, given no further info, only knows 1st dom
    end
    
  end % methods
end
