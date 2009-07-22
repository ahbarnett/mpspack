% CCAVITYMETALSCATTDEMO - Exponentially accurate computation of scattering
% from resonant polygonal C-shaped cavity.  Barnett 7/22/09.  Attempt 2

clear all
k = 3;%9.55;      % Wavenumber (around k=8.9 resonant for w=.12, ic=.2+.3i)
% geom of scatterer: ...based on unit square (-1/2,1/2)^2
ic = 0.2 + 0.3i; % inside corner
w = 0.08;%12;         % half-width of opening
% subdomain (mesh) params:
r = 0.95;               % Radius of outer circle
rmfs = 1/sqrt(2);      % MFS circle
c = 0.24; e = 0.3;     % half-size of subdomain around exit corner (c=horiz)
M = 40; %ceil(6*k);   % Number of quadrature points on O(1)-length segments
N = 100; %ceil(15*k);   % overall basis size scaling (M=30, N=70 for k=3)

% set up all segments:
icp = -.5+c + 1i*(w + (imag(ic)-w)*c/(imag(ic)+.5)); % internal corner pt
s = segment.polyseglist(M, [.5, .5+.5i, .5i, -.5+.5i, -.5+1i*e, -.5+1i*w, ...
                    icp, ic, real(ic)]);  % upper half of C-shaped scatt
s([5 6]).requadrature(M/2); % # quad pts on `short' segs near exit
s = [s(1:end-1) segment(M, [.5i 1i*r]) segment(M/2, [-.5+1i*e -.5-c+1i*e]) ...
     segment(M/2, [icp real(icp)]) ];   % 3 artificial bdries
al = asin(e/r);  % angle where artificial horiz line hits r-circle
s = [s segment(2*M, [0 r 0 pi/2]) segment(2*M, [0 r pi/2 pi-al]) ...
    segment(M/2, [0 r pi-al pi]) ]; % 3 arcs
s = [s segment(M, [-.5-c+1i*e r*exp(1i*(pi-al))]) ...
     segment(M, [-.5-c+1i*e -.5-c]) ];

s = [s s.reflect segment(M, [.5 r]) segment(M, [-.5+c real(ic)]) ...
     segment(M, [-.5-c -.5+c])]; % add reflections and x-axis segs
sdecomp=s([9:16 25:numel(s)]);   % All artificial boundaries (Gamma_ij)
sq=s([1:8 16+(1:8)]);            % All segments that are part of the square
sext=s([12:14 30 29 28]);           % ext bdry, Gamma_ie in CCW linking order
figure; s.plot; figure; sdecomp.plot;  % plot three types of segments in colors
h=sq.plot; set(h,'color', 'g'); h=sext.plot;set(h,'color', 'r'); 

ext = domain([], [], sext(end:-1:1), [1 1 1 -1 -1 -1]); % the exterior domain
d(1) = domain(s([33 12 9 2 1]),[1 1 -1 -1 -1]);    % 8 subdomains in CCW sense
d(2) = domain(s([9 13 15 10 4 3]),[1 1 -1 -1 -1 -1]);
d(3) = domain(s([16 35 11 6 5 10]),[1 1 -1 -1 -1 1]);
d(4) = domain(s([34 8 7 11]),[1 -1 -1 1]);
d(5) = domain(s([33 17 18 25 28]),[-1 1 1 1 -1]);
d(6) = domain(s([25 19 20 26 31 29]),[-1 1 1 1 1 -1]);
d(7) = domain(s([32 26 21 22 27 35]),[-1 -1 1 1 1 -1]);
d(8) = domain(s([34 27 23 24]),[-1 -1 1 1]);
dr = domain(s([15 14 30 31 32 16]),[1 1 -1 -1 1 -1]);  % regular domain
%figure; d.plot;   % shows each corner domain a different color

sdecomp.setmatch([k -k],[1 -1]); % Matching conditions for artificial boundaries

Z=@(t) rmfs*exp(2i*pi*t); Zp=@(t) 2i*pi*Z(t);   % MFS basis
ext.addmfsbasis({Z, Zp}, N, struct('eta',k, 'fast',2, 'nmultiplier',1.0));
o.type='s';                      % set up sine-type for each corner (sound-soft)
whichcorn = [5 6 5 3 3 3 4 4];   % which corner of each domain gets a nufbbasis?
for i=1:numel(d)
  j = whichcorn(i); cang = d(i).cang;  % choose N via corner angle and
  ra = max(abs(d(i).cloc(j)-d(i).x)); o.rescale_rad = 0.9*ra; % subdomain radius
  o.cornermultipliers = ((1:numel(cang))==j)*cang(j)/pi .* ra;
  d(i).addcornerbases(N, o);
end
o.nmultiplier = .2; o.rescale_rad = e; dr.addregfbbasis(-.75-c/2, N, o); % reg

pr=scattering(ext, [d dr]); %figure; pr.plot; % setup problem & show everything!
pr.setoverallwavenumber(k); pr.setincidentwave(-pi/6);
pr.updateN(N); diff([pr.basnoff pr.N])          % spits out bas{:}.Nf

% Solve and print solution
tic; pr.solvecoeffs; fprintf('\tcoeffs done in %.2g sec\n', toc)
fprintf('\tL2 bdry error norm = %g, coeff norm = %g\n', ...
        pr.bcresidualnorm, norm(pr.co))
figure; plot(real(pr.A*pr.co - pr.rhs)); title('residual vector');
if 1, o.bb=1.3*[-1 1 -1 1]; o.dx=0.03;
tic; [ui gx gy] = pr.gridincidentwave(o); u = pr.gridsolution(o); toc;
figure; imagesc(gx, gy, real(ui+u)); c = caxis; caxis([-1 1]*max(c));
axis equal tight; colorbar; set(gca,'ydir','normal'); end

