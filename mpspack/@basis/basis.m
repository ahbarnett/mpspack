classdef basis < handle
    
    % Class basis - Abstract class that defines the interfaces which are
    % common for all basis objects.
    
    properties
        k % Wave vector
        N % Number or degree of basis fct. Exact specification depends on
          % the type of basis
    end 
    methods (Abstract)
        [A, A1, A2] = eval(b,pts) 
                                % Evaluate a basis on a set of points
                                % A=eval(pts) returns only fct. values
                                % [A An]=eval(pts) returns fct. values plus
                                % normal derivatives.
                                % [A Ax Ay] returns fct. values plus x,y
                                % derivatives
                                
    end
end
                                
    
    