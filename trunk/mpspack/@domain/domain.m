% DOMAIN - create an interior domain possibly with excluded subregions
%
% d = DOMAIN(s, pm) creates an interior domain whose boundary is the
% segment array s
%
%
% d = DOMAIN(s, pm, si, pmi) si and sensi are cell arrays
%
% d = DOMAIN(s, pm, [], [], k)
%
% input a basis too.

classdef domain < handle
    properties
        segs                      % pointers to segments (row vec)
        pm                        % sense (+-1) of each segment (row vec)
%        csegs                    % numel(cor)-by-2 corner seg connections
        angoffs                   % Offsets for corner angles (complex #)
        angs                      % corner angles facing domain (complex #)
        cloc                      % corner locations (NaN if not a corner)
        k                         % wavevector for domain
        % n                         % refractive index of domain
    end
    methods
      function d = domain(s, pm, k)
        if nargin<3, k=NaN; end              % default wavenumber
        d.segs = s; d.pm = pm;
        i = (1-pm)/2+1; % indices (1,2) of which end starts the segment
                        % ahead of each corner
        % now make versions of the above 3 lists for the previous segment
        prevs = circshift(s,[0 1]); prevpm = circshift(pm,[0 1]);
        previ = (1+prevpm)/2+1; % ind (1,2) for segment behind each corner
        eps = 1e-15;                % max allowed corner position error
        for j=1:length(s)
          d.cloc(j) = s(j).eloc(i(j));
          if abs(d.cloc(j) - prevs(j).eloc(previ(j))) > eps
            disp(sprintf('domain warning: corner %d not connected!', j));
            d.cloc(j) = NaN;
          end
          d.angoffs(j) = s(j).eang(i(j));
          d.angs(j)    = -prevs(j).eang(previ(j))/d.angoffs(j);
        end
      end

        function h = plot(d)
          h = domain.showsegments(d.segs, d.pm);
          hold on; l=0.1;
          for j=1:numel(d.cloc)
            plot(real(d.cloc(j)), imag(d.cloc(j)), '.k', ...
                 'markersize', 20);
            % patch to illustrate corner angle...
            
          end
        end
        
        function i = inside(d,p)
          v = domain.approxpolygon(d.segs, d.pm);
          i = inpolygon(real(p), imag(p), real(v), imag(v));
        end
        
        e = isexterior(d)
        deletecorner(d)
    end
    
    methods(Static)    % invisible to outside, don't req domain
                       % obj to exist to call them
      v = approxpolygon(segs, pm)
      h = showsegments(segs, pm)
    end
end
