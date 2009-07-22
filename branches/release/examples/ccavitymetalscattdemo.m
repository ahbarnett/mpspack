% CCAVITYMETALSCATTDEMO - Exponentially accurate computation of scattering
% from resonant polygonal C-shaped cavity.  Barnett 7/22/09

clear all
k = 8.8;      % Wavenumber (around k=8.8 resonant for w=.12, ic=.2+.3i)
% geom of scatterer: ...based on unit square (-1/2,1/2)^2
ic = 0.2 + 0.3i; % inside corner
w = 0.12;         % half-width of opening
% subdomain (mesh) params:
r = 1.0;               % Radius of outer circle
rmfs = 1/sqrt(2);      % MFS circle
c = 0.15; e = 0.2;     % half-size of subdomain around exit corner (c=horiz)
M = 30; % Number of quadrature points on O(1)-length segments
N = 5*k; % overall basis size scaling

% set up all segments:
icp = -.5+c + 1i*(w + (imag(ic)-w)*c/(imag(ic)+.5)); % internal corner pt
s = segment.polyseglist(M, [.5, .5+.5i, .5i, -.5+.5i, -.5+1i*e, -.5+1i*w, ...
                    icp, ic, real(ic)]);  % upper half of C-shaped scatt
s([5 6]).requadrature(M/2); % # quad pts on `short' segs near exit
s = [s(1:end-1) segment(M, [.5i 1i*r]) segment(M/2, [-.5+1i*e -.5-c]) ...
     segment(M/2, [icp real(icp)]) ];   % 3 artificial bdries
s = [s segment(2*M, [0 r 0 pi/2]) segment(2*M, [0 r pi/2 pi]) ]; % 2 arcs

s = [s s.reflect segment(M, [.5 r]) segment(M, [-.5+c real(ic)]) ...
     segment(M, [-.5-c -.5+c]) segment(M, [-r -.5-c]) ]; % add x-axis segs
sdecomp=s([9:13 22:numel(s)]);   % All artificial boundaries (Gamma_ij)
sq=s([1:8 13+(1:8)]);            % All segments that are part of the square
sext=s([12 13 26 25]);           % ext bdry, Gamma_ie in CCW linking order
%figure; s.plot; figure; sdecomp.plot;  % plot three types of segments in colors
%h=sq.plot; set(h,'color', 'g'); h=sext.plot;set(h,'color', 'r'); 

ext = domain([], [], sext(end:-1:1), [1 1 -1 -1]); % the exterior domain
d(1) = domain(s([27 12 9 2 1]),[1 1 -1 -1 -1]);    % 8 subdomains in CCW sense
d(2) = domain(s([9 13 30 10 4 3]),[1 1 1 -1 -1 -1]);
d(3) = domain(s([29 11 6 5 10]),[1 -1 -1 -1 1]);
d(4) = domain(s([28 8 7 11]),[1 -1 -1 1]);
d(5) = domain(s([27 14 15 22 25]),[-1 1 1 1 -1]);
d(6) = domain(s([22 16 17 23 30 26]),[-1 1 1 1 -1 -1]);
d(7) = domain(s([29 23 18 19 24]),[-1 -1 1 1 1]);
d(8) = domain(s([28 24 20 21]),[-1 -1 1 1]);
%figure; d.plot;   % shows each domain a different color

sdecomp.setmatch([k -k],[1 -1]); % Matching conditions for artificial boundaries

Z=@(t) rmfs*exp(2i*pi*t); Zp=@(t) 2i*pi*Z(t);   % MFS basis
ext.addmfsbasis({Z, Zp}, N, struct('eta',k, 'fast',2, 'nmultiplier',2.0));
o.type='s';                      % set up sine-type for each corner (sound-soft)
whichcorn = [5 6 4 3 3 3 4 4];   % which corner of each domain gets a nufbbasis?
for i=1:numel(d)
  j = whichcorn(i); cang = d(i).cang;  % choose N via corner angle and
  ra = max(abs(d(i).cloc(j)-d(i).x)); o.rescale_rad = 0.9*ra; % subdomain radius
  o.cornermultipliers = ((1:numel(cang))==j)*cang(j)/pi .* ra;
  d(i).addcornerbases(N, o);
end

pr=scattering(ext, d); %figure; pr.plot;   % setup problem & show everything!
pr.setoverallwavenumber(k); pr.setincidentwave(-pi/6);

% Solve and print solution
tic; pr.solvecoeffs; fprintf('\tcoeffs done in %.2g sec\n', toc)
fprintf('\tL2 bdry error norm = %g, coeff norm = %g\n', ...
        pr.bcresidualnorm, norm(pr.co))
figure; plot(real(pr.A*pr.co - pr.rhs)); title('residual vector');
if 0, o.bb=1.3*[-1 1 -1 1]; o.dx=0.03;
tic; [ui gx gy] = pr.gridincidentwave(o); u = pr.gridsolution(o); toc;
figure; imagesc(gx, gy, real(ui+u)); c = caxis; caxis([-1 1]*max(c));
axis equal tight; colorbar; set(gca,'ydir','normal'); end

