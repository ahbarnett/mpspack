% POINTFARFIELD - evaluate far field from all bases at specified points
%
% u = p.pointfarfield(pts) computes the farfield u at the points pts.
%
% Copyright (C) 2014 Stuart C. Hawkins

function u = pointfarfield(self,pts)

% Note: no blocking implemented... typically number of points is much
% smaller than for two dimensional plots of the solution because the far
% field is computed in one dimension

% initialise the field
u = zeros(length(pts),1);

% set up points object for evalfarfield
points = struct('x',pts);

% loop through the domains
for n=1:numel(self.doms)
        
    % check if the domain is an exterior domain... non-exterior domains do
    % not contribute to the far field
    if self.doms(n).exterior
        
        % get hold of the domain
        dom = self.doms(n);
        
        % loop through bases
        for i=1:numel(self.bas)
            
            % get hold of the basis
            basis = self.bas{i};

            % check that the basis corresponds to the domain
            if utils.isin(basis,dom.bas);
            
                % evaluate the far field from the basis
                A = basis.evalfarfield(points);
                
                % get the coefficients
                cof = self.co(self.basnoff(i)+(1:basis.Nf),:);
                
                % add contribution to the far field
                u = u + A*cof;
                
            end
            
        end % end loop through bases
        
    end       
    
end % end loop through domains