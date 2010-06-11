% SCATTERING - define a frequency-domain scattering problem
%
%  pr = scattering(airdoms, doms) creates a scattering problem
%   object pr which comprises the `air domains' (exterior parts of domain which
%   will have incident field in them) airdoms, and other interior domains doms
%   which will have zero incident field. airdoms and doms are row vecs of
%   domain handles, although doms may be empty (scattering off metallic object).
%
% See also BVP, PROBLEM

classdef scattering < bvp & handle
  properties
    incang                             % incident angle
    ui                                 % incident field function in air domain
    uix, uiy                           % x,y-deriv functions of above
  end
  
  methods % ---------------------------------- methods particular to scattering
    function pr = scattering(airdoms, doms) % .............. constructor
      if nargin==0, airdoms = []; doms = []; end     % needs empty constructor
      pr = pr@bvp([airdoms doms]);      % make an instance of a BVP (see Docu)
      if ~isempty([airdoms doms])       % if not an empty constructor...
        if numel(airdoms)==0, error('must be at least one air domain'); end
        for d=airdoms, d.isair = 1; end   % flag the domains as air or not
        for d=doms, d.isair = 0; end
        if max(vertcat(airdoms.refr_ind))~=min(vertcat(airdoms.refr_ind))
          fprintf('warning: all air domains need same refractive index!\n')
        end
      end
    end
  
    function setincidentwave(pr, t)
    % SETINCIDENTWAVE - choose plane wave or u_incident field, compute f,g BCs
    %
    %  setincidentwave(pr, t) sets up the inhomogeneous BCs required in
    %   scattering problem pr, for planewave at angle t in [0,2pi].
    %   Wavenumber must already be chosen by the problem k.
    %
    %  setincidentwave(pr, ui, uix, uiy) uses a given incident field function
    %   ui and its x,y-deriv functions uix, uiy. The functions take a complex
    %   number z = x + iy where (x,y) is a point in the plane. Warning: in
    %   this case it is up to the user to ensure ui is a Helmholtz solution of
    %   the correct wavenumber in the air domain(s)!
    %
    %  Careful: calling this routine overwrites all inhomogeneity functions or
    %   data f, g stored on any of the problem's segments. However it preserves
    %   existing a, b BC or matching coeffs (which must be set up on entry).
      if nargin==2
        if isempty(pr.k), error('please set wavenumber before choosing incident plane wave angle!'); end
        kvec = pr.k*exp(1i*t);                % set up a plane-wave field
        pr.incang = t;
        ui = @(x) exp(1i*real(conj(kvec) .* x));
        uix = @(x) 1i*real(kvec)*ui(x); uiy = @(x) 1i*imag(kvec)*ui(x);
      end
      pr.ui = ui; pr.uix = uix; pr.uiy = uiy;
      % now set up inhomogeneities in BCs, incident air field & no other field..
      for s=pr.segs
        if s.bcside==0                       % matching, may be dielectric
        if s.dom{1}.isair==s.dom{2}.isair    % either both incident or both zero
          s.f = @(t) zeros(size(s.t)); s.g = s.f;    % homogeneous matching
        else
          % extract the relevant coeffs for whichever side have inc field on it:
          if s.dom{1}.isair       % + side has inc, - side has zero field
            a = s.a(1); b = s.b(1);
          else                    % - side has inc, + side has zero field
            a = s.a(2); b = s.b(2);
          end
          s.f = @(t) -a * ui(s.Z(t));        % jump in value of u_s (NB minus)
          s.g = @(t) -b * (uix(s.Z(t)).*real(s.Zn(t)) + ...
                uiy(s.Z(t)).*imag(s.Zn(t))); % jump in deriv of u_s
        end
        elseif s.bcside==1 | s.bcside==-1    % BC
          ind = (1-s.bcside)/2+1; % index 1 or 2 for which side the BC on
          d = s.dom{ind};         % handle of domain on the revelant side
          if d.isair              % air-to-metallic boundary
            if s.b==0             % Dirichlet only. Minus sign: u_s cancels u_i
              s.f = @(t) -s.a * ui(s.Z(t));  % NB keeping f a func not an array
            else                  % Robin, Neumann
              s.f = @(t) -s.a * ui(s.Z(t)) -s.b * (uix(s.Z(t)).*real(s.Zn(t))...
                              + uiy(s.Z(t)).*imag(s.Zn(t)));
            end
          else                    % internal-to-metallic boundary
            s.f = @(t) zeros(size(s.t));     % f = homogeneous BCs
          end
        end
      end % segs loop
    end % func
    
    function [u di] = pointincidentwave(pr,p,o) % ........eval u_inc on pointset
    % POINTINCIDENTWAVE - evaluate incident wave for a problem on a pointset
    %
    %  [ui di] = pointincidentwave(pr, pts) returns array of values ui, and
    %   optionally, domain index list di (integer array of same shape as ui).
    %   Decisions about which domain a gridpoint is in are done using
    %   domain.inside, which may be approximate. di and ui are NaN outside of
    %   all domains.
    %
    % See also GRIDINCIDENTWAVE
    if nargin<3, o=[]; end
    if ~isfield(o,'all'), o.all=0; end; % Evaluate wave over all domains
      di = NaN*zeros(size(p.x));                    % NaN indicates in no domain
      u = di;                                       % solution field
      for n=1:numel(pr.doms)
        d = pr.doms(n);
        ii = d.inside(p.x);
        di(ii) = n;
        if d.isair | o.all                          % here there's u_i wave
          u(ii) = pr.ui(p.x(ii));
        else                                        % here there's 0 inc wave
          u(ii) = 0;
        end
      end
    end % func
    
    function [u gx gy di] = gridincidentwave(pr, o) % ...... eval u_inc on grid
    % GRIDINCIDENTWAVE - evaluate u_inc over a grid
    %
    %  [ui gx gy di] = gridincidentwave(pr, opts) returns array of values ui,
    %   optionally, x- and y-grids (1d lists) gx, gy, and domain index list di
    %   (integer array of same shape as ui). Decisions about which domain a
    %   gridpoint is in are done using domain.inside, which may be approximate.
    %
    % See also POINTINCIDENTWAVE, GRIDBOUNDINGBOX
    
      o = pr.gridboundingbox(o);
      gx = o.bb(1):o.dx:o.bb(2); gy = o.bb(3):o.dx:o.bb(4);  % plotting region
      [xx yy] = meshgrid(gx, gy); zz = xx(:) + 1i*yy(:);  % keep zz rect array
      [u di] = pr.pointincidentwave(pointset(zz),o);
      u=reshape(u,size(xx,1),size(xx,2));
      di=reshape(di,size(xx,1),size(xx,2));      
    end % func
    
    function showthreefields(pr, o)
    % SHOWTHREEFIELDS - plot figure with u_i, u_s and u_t on subplot grids (Re)
    %
    %   if opts.testtransparent is true, then the final plot is error from
    %    transparency, reports L2 error estimated over grid too.
    %   if opts.imag = true, plots imag instead of real part
    %   opts.bdry = true, shows boundary too
    %   opts.sepfigs: if true, make three separate figures
    %
      if nargin<2, o = []; end
      if ~isfield(o, 'imag'), o.imag = 0; end
      if ~isfield(o, 'bdry'), o.bdry = 0; end
      if ~isfield(o,'sepfigs'), o.sepfigs=0; end % If true plot in sep. figs.
      if o.sepfigs,
          figure;
      else
          tsubplot(1,3,1);
      end
      [ui gx gy di] = pr.gridincidentwave(o);
      if o.imag, imagesc(gx, gy, imag(ui)); title('Im[u_i]');
      else, imagesc(gx, gy, real(ui)); title('Re[u_i]'); end
      c = caxis; caxis([-1 1]*max(c));           % make colorscale symmetric
      axis equal tight;colorbar; set(gca,'ydir','normal'); hold on;
      if o.bdry, pr.showbdry; end
      
      if o.sepfigs,
          figure
      else
          tsubplot(1,3,2);
      end
      [us gx gy di] = pr.gridsolution(o);
      if o.imag, imagesc(gx, gy, imag(us)); title('Im[u_s]');
      else, imagesc(gx, gy, real(us)); title('Re[u_s]'); end
      %c = caxis; caxis([-1 1]*max(c));
      caxis(2*[-1 1]*max(c));                 % choose double caxis hack
      axis equal tight;colorbar; set(gca,'ydir','normal'); hold on;
      if o.bdry, pr.showbdry; end
      
      if o.sepfigs,
          figure
      else
        tsubplot(1,3,3);
      end
      if isfield(o,'testtransparent') & o.testtransparent
        [xx yy] = meshgrid(gx, gy); zz = xx + 1i*yy;
        uiR = pr.ui(zz);                  % u_i eval over whole R^2
        uerr = ui+us-uiR;
        if o.imag
          imagesc(gx, gy, imag(uerr));title('Im[u_{err}] = Im[u - u_i(R^2)]');
        else
          imagesc(gx, gy, real(uerr));title('Re[u_{err}] = Re[u - u_i(R^2)]');
        end
        fprintf('L2 transparency error est on grid = %g\n', o.dx*norm(uerr))
      else
        if o.imag, imagesc(gx, gy, imag(ui+us)); title('Im[u] = Im[u_i+u_s]');
        else, imagesc(gx, gy, real(ui+us)); title('Re[u] = Re[u_i+u_s]'); end
      end
      %c = caxis; caxis([-1 1]*max(c));
      caxis(2*[-1 1]*max(c));                  % choose double caxis hack
      axis equal tight;colorbar; set(gca,'ydir','normal'); hold on;
      if o.bdry, pr.showbdry; end
    end % func
    
    function showfullfield(pr, o)
    % SHOWFULLFIELD - plot u_t = u_i+u on the grid
    %
    %   opts.imag = true, plots imag instead of real part
    %   opts.bdry = true, shows boundary too

    % Notes: Timo's idea.
      if nargin<2, o = []; end
      if ~isfield(o, 'imag'), o.imag = 0; end
      if ~isfield(o, 'bdry'), o.bdry = 0; end
    
      [ui gx gy di] = pr.gridincidentwave(o);
      u=pr.gridsolution(o);
      u=ui + u;
      figure;
      if o.imag
          imagesc(gx, gy, imag(u));title('Im[u_{tot}]');
      else
          imagesc(gx, gy, real(u)); title('Re[u_{tot}]');
      end
      c = caxis; caxis([-1 1]*max(c));
      axis equal tight;colorbar; set(gca,'ydir','normal'); hold on;
      if o.bdry, pr.showbdry; end
    end      
      
  end % methods
end
