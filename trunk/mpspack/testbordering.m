% test segment-domain bordering
% 7/10/08

clear all classes
verb = 1;
% one analytic segment, periodic radial function
a = 0.3; w = 3;
M = 100;
R = @(t) 1 + a*cos(w*t); Rt = @(t) -w*a*sin(w*t);
Z = @(s) exp(2i*pi*s).*R(2*pi*s);
sa = segment(M, {Z, @(s) 2*pi*(1i*Z(s) + exp(2i*pi*s).*Rt(2*pi*s))}, 'p');
sq = segment.polyseglist(10, -1.7-0.5i + 0.5*[0 1i 1+1i 1]);  % CW square
sc = segment(30, [1+1.5i, 0.5, 0, 2*pi], 'p'); % CCW circle
d = domain([], [], {sa sq sc}, {-1 1 -1});
% connect up some interior domains with the above exterior domain
da = domain(sa, 1); dq = domain(sq(end:-1:1), -1); dc = domain(sc, 1);
if verb, figure; opts.gridinside=0.1; h=domain.showdomains([d da dq dc], opts);
   axis equal; title('test domain connections: ext w/ holes & 3 interior'); end
s = [sa sq sc];
disp('There should have been no errors, and the following should show two domain objects in each cell array:')
s.dom               % list each segment - how identify/distinguish domains ?
                    % by name? give symbol to each?
disp('The plot should show 3 interior and 1 exterior domain, of different colors. There should be one color of normal vector showing with each sense on eac boundary')
