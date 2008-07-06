classdef basis < handle
    
    % Class basis - Abstract class that defines the interfaces which are
    % common for all basis objects.
    
    properties
        k % Wave vector
        N % Number or degree of basis fct. Exact specification depends on
          % the type of basis
    end 
    methods (Abstract)
        [A, An, Ax, Ay] = eval(b,pts) 
                                % Evaluate a basis on a set of points
                                % A=eval(pts) returns only fct. values
                                % [A An Ax Ay] returns fct. values plus x,y
                                % derivatives
                                
    end
end
                                
    
    