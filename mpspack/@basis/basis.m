classdef basis < handle
    
    % Class basis - Abstract class that defines the interfaces which are
    % common for all basis objects.
    
    properties
        k                       % Wavevector
        N                       % Number or degree of basis fct. Exact
                                % specification depends on type of basis
        Nf                      % Actual number of basis functions
    end 
    methods (Abstract)
        [A, A1, A2] = eval(b,pts) 
                                % Evaluate a basis on a set of points
                                % A=eval(pts) returns only fct. values
                                % [A A1]=eval(pts) returns fct. values plus normal
                                % derivatives.
                                % [A A1 A2]=eval(pts) returns fct. values plus
                                % x and y derivatives.
                                
    end
end
                                
    
    