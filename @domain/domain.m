% DOMAIN - create an interior/exterior domain possibly with excluded subregions
%
% A domain is an ordered connected list of segments defining the exterior
% boundary, with zero or more ordered connected lists of segments defining the
% boundaries of any interior excluded regions. If the exterior boundary is empty
% the domain from which interior regions are possibly excluded is taken to be
% the whole plane, resulting in an unbounded exterior domain.
%
%  d = DOMAIN(s, pm) creates an interior domain whose boundary is the list
%   of handles of segment objects s, using the list of senses pm (each element
%   is the number +1 or -1). A warning is given if the segments do not appear
%   to connect up at corners. Normals should all point outwards, ie away from
%   the domain, otherwise a warning is given. If pm has only 1 element, it will
%   be duplicated as necessary to match the size of s.
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
% See also: SEGMENT, domain/PLOT

% Copyright (C) 2008 - 2012, Alex Barnett, Timo Betcke


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
        exterior                  % true if exterior domain (w bndd complement) 
        bas                       % pointer list to basis set objects (cell arr)
        k                         % wavevector for domain: controls all bases
        refr_ind                  % refractive index of domain (default = 1)
        isair                     % 1 if gets inc wave; 0,2,3.. if not (scatt)
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
        d.refr_ind = 1.0;                       % default refractive index
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
           p(j) = d.seg(j).Z(t) + eps * d.pm(j) * d.seg(j).Zn(t);
        end
        norout = ~d.inside(p);                 % true if points away from domain
      end
      
      function i = inside(d, p) % .................................... inside
      % INSIDE - return true (false) for points inside (outside) a domain
      %
      % i = INSIDE(d, p) returns a logical array of values specifying if
      %   each point in the list of complex numbers p is inside (1) or outside
      %   (0) the domain d. The output i is the same shape as input p.
      %
      % Issues/notes:
      %  * Uses fast inpolygon implementations in inpolywrapper
      %  * inluded temporary hack to make semi-infinite strip domain if corners
      %    don't connect up - the domain is officially borken then anyway!
        if d.exterior
          i = logical(ones(size(p(:))));
        else                                % interior domain
          js = find(d.spiece==0);           % indices of segs outer bdry piece
          v = domain.approxpolygon(d.seg(js), d.pm(js));
          i = utils.inpolywrapper(p(:), v);
        end
        for piece=1:max(d.spiece)         % kill pts from each interior piece
          js = find(d.spiece==piece);
          v = domain.approxpolygon(d.seg(js), d.pm(js));
          % HACK which extends an excluded strip down/upwards:  
          if isnan(d.cloc), e=(d.pm(js)+3)/2; v = [d.seg(js).eloc(3-e)+d.pm(js)*10i; v; d.seg(js).eloc(e); d.seg(js).eloc(e)+d.pm(js)*10i]; end
          i = i & ~utils.inpolywrapper(p(:), v);
        end
        i = reshape(i, size(p));
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
      
      function sp = speed(d) % ........... get all quadr speeds assoc w/ domain
      % SPEED - return col vector of quadrature speeds for a domain boundary
        sp = [];
        for j=1:numel(d.seg)
          if d.pm(j)==1
            sp = [sp; d.seg(j).speed];
          else                            % reverse order
            sp = [sp; d.seg(j).speed(end:-1:1)];
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
      % BOUNDINGBOX - return bounding [xmin xmax ymin ymax] for interior domain
      %
      %  If the domain is exterior, a box enclosing the boundary with some
      %   extra padding is returned (padding currently not user-selectable).
        x = d.x;
        bb = [min(real(x)), max(real(x)), min(imag(x)), max(imag(x))];
        if isempty(bb), bb = zeros(1,4); end    % only if d is whole plane
        if isnan(d.exterior) || d.exterior
          pad = 0.5;                      % pad exterior region for ext domains
          bb = bb + pad*[-1 1 -1 1];
        end
      end
      
      function xc = center(d) % .............. center of bounding box
      % CENTER - return center x (as complex number) of bounding box of domain
      %
      % c = CENTER(dom) returns center of domain dom's rectangular bounding box.
      %
      % See also: DOMAIN.BOUNDINGBOX
        bb = boundingbox(d);
        xc = (bb(1)+bb(2)+1i*bb(3)+1i*bb(4))/2;
      end
      
      function diam = diam(d) % ..... diameter of interior domain (about center)
      % DIAM - approximate radius of an interior domain
      %
      %  d = diam(dom) returns the radius of domain dom, defined as the
      %   maximum distance of any of its defining points from its `center'.
      %
      % See also: DOMAIN.CENTER
        diam = max(abs(d.x - center(d))); % not quite optimal since uses box
      end
      
      function [zz ii gx gy] = grid(d, dx, bb) % ......... grid covering domain
      % GRID - make grid covering interior domain, or some of exterior domain
      %
      %  [zz ii gx gy] = grid(dom, dx) returns column-vector list of C-#s zz
      %   falling in the domain, and the indices ii to where these fall in a
      %   regular rectangular grid with x- an y-axis 1D grids gx and gy
      %   respectively.
      %
      %  [zz ...] = grid(dom, dx, bb) overrides bounding box
      %
      %  See also: DOMAIN.BOUNDINGBOX
        if dx<=0, error('dx must be positive!'); end
        if nargin<3, bb = d.boundingbox; end
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
        for i=1:numel(d.bas); Nf = Nf + d.bas{i}.Nf; end
      end
      
      function clearbases(d) % .............. removes all basis sets from domain
      % CLEARBASES - remove all basis set associations from a domain
        d.bas = {};    % possibly we should also kill any d.bas{:}.doms == d
      end
      
      function showbasesgeom(d) % ................. show geometry of basis objs
        for i=1:numel(d.bas)
          opts.label = sprintf('%d', i);              % label by domain's bas #
          d.bas{i}.showgeom(opts);
        end
      end
      
      function [c h mis gx gy] = showimagparam(d, o)  % ..... invert Z imag dist
      % SHOWIMAGPARAM - compute and show imag param dist on grid in domain
      %
      % [c h mis gx gy] = showimagparam(d, o) adds contour lines of constant
      %   Im part of S(z) computed using z on a grid filling the interior of
      %   domain d.  mis returns an array of the imag part (or nan if outside),
      %   on the product grid given by gx and gy.
      %
      % o.dx : grid spacing (smaller dx is slower)
      % o.bb : override bounding box
      % o.levels : levels to feed to contouring
      % o.outside : 0 (default) shows only inside; 1 shows inside & outside.
      %
      % Other options fed to segment.invertZparam
        if nargin<2, o = []; end
        if ~isfield(o, 'dx'), o.dx = 0.03; end       % default opts
        if ~isfield(o, 'outside'), o.outside = 0; end  % default: inside only
        if isfield(o, 'bb'); [zz ii gx gy] = d.grid(o.dx, o.bb);
        else, [zz ii gx gy] = d.grid(o.dx); end
        if o.outside, [xx yy] = meshgrid(gx,gy);
          zz = xx(:) + 1i*yy(:); ii = true(size(xx)); end % plot inside&outside
        s = d.seg;
        if numel(s)~=1, error('domain must have exactly 1 segment!'); end
        t = s.invertZparam(zz,o);      % expensive part
        if ~o.outside, t(find(imag(t)<0)) = nan+1i*nan; end % note imag(nan)=0 so need 1i*nan too!
        mis = nan*ii;       % output array
        mis(ii) = min(abs(imag(t)), [], 1);
        if isfield(o, 'levels'), L = o.levels; else, L = 0:0.01:0.1; end
        [c h] = contour(gx, gy, mis, L, 'k-', 'linewidth', 1.5); colorbar;
      end
      
      % methods defined by separate files...
      addconnectedsegs(d, s, pm, o)      % helper routine for constructor
      h = plot(d, o)                     % domain plot: o is plot opts struct
      addregfbbasis(d, varargin)         % add reg FB basis object
      addnufbbasis(d,varargin)           % add irreg. FB basis
      addcornerbases(d, N, opts)         % add multiple nu-FB's at corners
      addrpwbasis(d, varargin)           % add real PW basis
      addmfsbasis(d, varargin)           % add MFS basis
      addqprayleighbasis(d, seg, pm, varargin) % add QP Rayleigh basis
      b = addlayerpot(d, segs, a, opts)  % add layer-potential basis
      [A A1 A2] = evalbases(d, p, opts)  % evaluate all basis funcs in domain
      setrefractiveindex(doms, n)
      
      % ****** not yet implemented ***** ( low priority, mostly)
      checktopology(d)             % checks all pieces in interior, normals
      requadrature(d, M, qtype)    % resets M on all segments, quadr type
      subtract(d, dint)            % remove (exclude) domain dint from d
    end % methods
    
    % --------------------------------------------------------------------
    methods(Static)    % these don't need domain obj to exist to call them...
      v = approxpolygon(seg, pm)
      h = showsegments(seg, pm, o)
      [x nx] = stackquadpts(seg, pm)
      h = showdomains(dlist, opts)
    end
end
