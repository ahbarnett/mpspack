% DOMAIN - create an interior/exterior domain possibly with excluded subregions
%
% A domain is an ordered connected list of segments defining the exterior
% boundary, with one or more ordered connected lists of segments defining the
% boundaries of excluded regions. If the exterior boundary is empty it is
% taken to be the whole plane, ie an exterior domain. We implement it as a
% handle class.
%
%  d = DOMAIN(s, pm) creates an interior domain whose boundary is the list
%   of handles of segment objects s, using the list of senses pm (each element
%   is the number +1 or -1). A warning is given if the segments do not appear
%   to connect up at corners. Normals should all point outwards, ie away from
%   the domain, otherwise a warning is given. If pm has only 1 element, it will
%   be duplicated as necessary.
%
%  d = DOMAIN() or d = DOMAIN([], []) creates an exterior domain equal to the
%   whole plane R^2.
%
%  d = DOMAIN([], [], si, pmi) creates an exterior domain equal to the whole
%   plane minus an excluded region whose (non-intersecting) boundary is given
%   by segment list si and sign list pmi (which have the same format as s, pm
%   above). If si and pmi are instead cell arrays of segment lists and
%   corresponding sign lists, each cell element is an excluded region. As
%   above, warnings are given if intersections or incorrect normals are found.
%
%  d = DOMAIN(s, pm, si, pmi) combines the above features, creating a bounded
%   domain with excluded region(s).
%
% *** Not yet implemented :
%  d = DOMAIN(s, pm, si, pmi, k) makes a domain with wavenumber k.
%
% *** Issues:
% * should diam, center, x, w, boundingbox, etc, be precomputed on construction,
%   only recomputed if a segment changes? Currently not. Not a big deal.
%
% See also: SEGMENT, domain/PLOT

classdef domain < handle
    properties
        seg                       % pointer list to segments (row vec)
        pm                        % sense (+-1) of each segment (row vec)
        spiece                    % which conn. piece each segment in (row vec)
        cangoff                   % offsets for corner angles (c-# row vec)
        cang                      % corner angs (0,2pi) facing domain (row vec)
        cloc                      % corner locations (NaN if not a corner)
        perim                     % perimeter (via quadr pts; spectral acc)
        area                      % area (via quadr pts; only O(1/M^2) acc)
        exterior                  % true if exterior domain (bnded complement) 
        bas                       % pointer list to basis set objects (cell arr)
        k                         % wavevector for domain
        noff                      % dof offset used as temp in bvp.fillbcmatrix
%       cseg                      % numel(cloc)-by-2 corner-seg connectivity
%       n                         % refractive index of domain
    end
    methods % ---------------------------------------------------------------

      function d = domain(s, pm, si, pmi, k, b) % .............. constructor
        if nargin<5, d.k=NaN; else, d.k = k; end      % set default wavenumber
        d.seg = []; d.pm = []; d.perim = 0;
        d.cloc = []; d.cang = []; d.cangoff = [];     % start corner lists
        if nargin>1 & ~isempty(s)              % seg list gives primary bdry
          d.area = 0; d.exterior = 0;
          addconnectedsegs(d, s, pm);
        else                                   % empty primary seg list
          d.area = Inf; d.exterior = 1; d.spiece = [];
        end
        if nargin>2 & ~isempty(si)             % parse interior segments if any
          if ~iscell(si)
            si = {si}; pmi = {pmi};            % if one int piece, make cells
          end
          for i=1:numel(si)                    % loop thru cell arr, int pieces
            dcur = utils.copy(d);              % deep copy of *current* domain
            oh = []; oh.hole = 1; addconnectedsegs(d, si{i}, pmi{i}, oh);
            js = find(d.spiece==i);            % topo test segments just added
            if ~isempty(find(~dcur.inside(vertcat(d.seg(js).x)))) % fell out?
              fprintf('domain warning: piece %d not in current domain!\n', i)
            end
          end
        end
        norout = normalscheck(d);              % normals sense of whole seg list
        if ~isempty(find(norout==0))           % any normals point wrong (in)?
          fprintf('domain warning: probable bad sense of segment normals!\n')
          norout
        end
        
        d.bas = {};                            % initialize w/ no basis sets
      end

      %function delete(d) % .............................. destructor
      % This removes links in the segment object connecting to this domain
      % However this is too messy. Use segment.disconnect
      %end
      
      function norout = normalscheck(d) % ....... check normal senses of domain
      % NORMALSCHECK - check that senses of normals point away from a domain
        eps = 1e-8;           % small dist from bdry to move in normal direc
        t = 1/2;              % where to test along seg; so Napprox must be even
        p = zeros(size(d.seg));
        for j = 1:numel(d.seg)
           p(j) = d.seg(j).Z(t) + eps * d.pm(j) * d.seg(j).normal(t);
        end
        norout = ~d.inside(p);                 % true if points away from domain
      end
      
      function i = inside(d, p) % .................................... inside
      % INSIDE - return true (false) for points inside (outside) a domain
        if d.exterior
          i = logical(ones(size(p)));
        else                                % interior domain
          js = find(d.spiece==0);           % indices of segs outer bdry piece
          v = domain.approxpolygon(d.seg(js), d.pm(js));
          i = inpolygon(real(p), imag(p), real(v), imag(v));
        end
        for piece=1:max(d.spiece)         % kill pts from each interior piece
          js = find(d.spiece==piece);
          v = domain.approxpolygon(d.seg(js), d.pm(js));
          i = i & ~inpolygon(real(p), imag(p), real(v), imag(v)); % NB logical
        end
      end

      function x = x(d) % ................. get all quadr pts assoc w/ domain
      % X - return column vector of quadrature points on a domain boundary
      x = domain.stackquadpts(d.seg, d.pm);
      end
      
      function nx = nx(d) % ............... get all normals assoc w/ domain
      % NX - return column vector of unit outward normals on a domain boundary
        nx = [];
        for j=1:numel(d.seg)
          if d.pm(j)==1
            nx = [nx; d.seg(j).nx];
          else                            % reverse order and negated
            nx = [nx; -d.seg(j).nx(end:-1:1)];
          end
        end
      end
      
      function w = w(d) % ................... get all quadr wei assoc w/ domain
      % W - return row vector of quadrature weights for a domain boundary
        w = [];
        for j=1:numel(d.seg)
          if d.pm(j)==1
            w = [w, d.seg(j).w];
          else                            % reverse order
            w = [w, d.seg(j).w(end:-1:1)];
          end
        end
      end
      
      function bb = boundingbox(d) % ......... bounding box
      % BOUNDINGBOX - return [xmin xmax ymin ymax] for interior domain
      %
      %  If the domain is exterior, a box enclosing the boundary with some
      %   extra padding is returned.
        x = d.x;
        bb = [min(real(x)), max(real(x)), min(imag(x)), max(imag(x))];
        if isempty(bb), bb = zeros(1,4); end    % only if d is whole plane
        if d.exterior
          pad = 0.5;                      % pad exterior region for ext domains
          bb = bb + pad*[-1 1 -1 1];
        end
      end
      
      function xc = center(d) % .............. center of bounding box
      % CENTER - return center x (as complex number) of bounding box of domain
        bb = boundingbox(d);
        xc = (bb(1)+bb(2)+1i*bb(3)+1i*bb(4))/2;
      end
      
      function diam = diam(d) % ..... diameter of interior domain (about center)
      % DIAM - approximate diameter of interior domain
        diam = max(abs(d.x - center(d))); % not quite optimal since uses box
      end
      
      function [zz ii gx gy] = grid(d, dx) % ............. grid covering domain
      % GRID - make grid covering interior domain, or some of exterior domain
      %
      %  [zz ii gx gy] = grid(dom, dx) returns list of C-#s zz and the indices
      %   ii to where these fall in a regular rectangular grid with x- and
      %   y-grids gx and gy respectively.
        if dx<=0, error('dx must be positive!'); end
        bb = d.boundingbox;
        bb(1) = dx * floor(bb(1)/dx);            % quantize to grid through 0
        bb(3) = dx * floor(bb(3)/dx);
        gx = bb(1):dx:bb(2); gy = bb(3):dx:bb(4);         % plotting region
        [xx yy] = meshgrid(gx, gy); zz = xx(:) + 1i*yy(:); ii = d.inside(zz);
        zz = zz(ii);                             % keep only inside points
        ii = reshape(ii, size(xx));              % return ii to rect array
      end
      
      function deletecorner(d,j)
      % DELETECORNER - remove a given corner number from a domain
        d.cloc(j) = NaN;                 % note this doesn't shorten c arrays
      end
      
      function Nf = Nf(d) % ................. total # of basis funcs in domain
      % NF - total number of basis functions (dofs) associated with domain
        Nf = 0;
        for b = d.bas, Nf = Nf + b{1}.Nf; end    % cell array {1} feels clumsy
      end
      
      function clearbases(d) % .............. removes all basis sets from domain
      % CLEARBASES - remove all basis sets from a domain
        d.bas = {};
      end
      
      % methods defined by separate files...
      addconnectedsegs(d, s, pm, o)  % helper routine for constructor
      h = plot(d, o)                     % domain plot: o is plot opts struct
      addregfbbasis(d, origin, N, k, opts) % add reg FB basis object
      addnufbbasis(d,origin,nu,offset,branch,N,k,opts) % add irreg. FB basis
      addrpwbasis(d, N, k, opts)           % add real PW basis
      addmfsbasis(d, Z, tau, N, k, opts)   % add MFS basis
      [A A1 A2] = evalbases(d, p)     % evaluate all basis funcs in domain
      
      % ****** not yet implemented ***** ( low priority, mostly)
      checktopology(d)             % checks all pieces in interior, normals
      requadrature(d, M, type, inds) % resets M, quadr, in some segments
      subtract(d, dint)              % remove (exclude) domain dint from d
      changek(d, k)                  % changes wave# of basis sets in domain
    end % methods
    
    % --------------------------------------------------------------------
    methods(Static)    % these don't need domain obj to exist to call them...
      v = approxpolygon(seg, pm)
      h = showsegments(seg, pm, o)
      x = stackquadpts(seg, pm)
      h = showdomains(dlist, opts)
    end
end
