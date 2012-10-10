% test the domain constructor, via plotting and computing area & perim.
% barnett 7/3/08. tweaked 10/8/12 to add 2 smooth segments making a hole.

clear all classes
verb = 1;                           % verbosity 0=quiet, 1=figures, 2=more

if verb, fprintf('You should see 9 figures with a variety of domains:\n'); end

% simply-connected polygon w/ different quadrature types
p = [0 1 1+2i]; M = 10;
s(1) = segment(M, [p(1) p(2)], 'c');
s(2) = segment(M, [p(2) p(3)], 't');
s(3) = segment(M, [p(1) p(3)], 't');
fprintf('triangle...   '), d = domain(s, [1 1 -1]);
if verb, figure; opts.gridinside=0.1; plot(d, opts); axis equal;
title('test domain: triangle');
fprintf('area err %g, perim err %g\n', d.area-1, d.perim-(3+sqrt(5))); end
if 0
  d = domain(s,ones(1,length(s)));   % should cause corner disconnect warning
end

% simply-connected lines + analytic curve
s(4) = segment(M, {@(s) (1-exp(1i*pi*s))/2, @(s) -1i*pi/2*exp(1i*pi*s)}, 'p');
fprintf('lines + analytic...\n'); disconnect(s);
d = domain([s(4); s(2); s(3)], [1 1 -1]);
if verb, figure; plot(d, opts); axis equal;title('test domain: lines + anal');
end

% simply-connected lines + analytic curve
s(4) = segment(M, [0.5 0.5 -pi 0]);
s(5) = segment(M, {@(t) 2i*(1-t) + 1 + 0.3*sin(2*pi*t), @(t) -2i+0.6*pi*cos(2*pi*t)});
fprintf('lines + arc + analytic...\n'); disconnect(s);
d = domain([s(4); s(5); s(3)], [1 -1 -1]);
if verb, figure; plot(d, opts); axis equal;title('test domain: mixed segs');
end

% one analytic segment, periodic radial function
a = 0.3; w = 3;
M = 100;
R = @(t) 1 + a*cos(w*t); Rt = @(t) -w*a*sin(w*t);
Z = @(s) exp(2i*pi*s).*R(2*pi*s);
sa = segment(M, {Z, @(s) 2*pi*(1i*Z(s) + exp(2i*pi*s).*Rt(2*pi*s))}, 'p');
fprintf('single radial func...\n')
d = domain(sa, [1]);
if verb, figure; plot(d, opts); axis equal; title('test domain: radial func');
end

% multiply-connected: genus 1
si(1) = segment(30, [.5 .3 0 2*pi], 'p');       % small CCW circle arc
fprintf('genus 1...\n'), disconnect(s);
d = domain([s(4); s(5); s(3)], [1 -1 -1], si, [-1]);  % place its reverse in d
if verb, figure; opts.gridinside=0.05; plot(d, opts); axis equal;
  title('test domain: multiply-connected genus 1'); end

% genus 2
p = [.5+0.6i, .9+1.6i, 1.2+1.4i];       % its area is 0.19
M = 10;
sii(1) = segment(M, [p(1) p(2)]);       % small CW tri 
sii(2) = segment(M, [p(2) p(3)]);
sii(3) = segment(M, [p(3) p(1)]);
fprintf('genus 2...\n'), disconnect(s); disconnect(si);
d = domain([s(4); s(5); s(3)], [1 -1 -1], {si, sii}, {[-1], [1 1 1]});
if verb, figure; h=plot(d, opts); axis equal; title('test domain: genus 2');end
if 0
  p(1) = [.3+0.1i];    % tri now intersects the interior circle: raises warning
  sii(1) = segment(M, [p(1) p(2)]);       % small CW tri 
  sii(3) = segment(M, [p(3) p(1)]);
  d = domain([s(4); s(5); s(3)], [1 -1 -1], {si, sii}, {[-1], [1 1 1]});
end

% whole plane
fprintf('whole plane...\n'), d = domain();
if verb, figure; opts.gridinside=0.1;
  h=plot(d, opts); axis equal; title('test domain: whole plane');end

% exterior domain
fprintf('exterior domain w/ 1 hole...\n'), disconnect(sa);
d = domain([], [], sa, -1);
if verb, figure; opts.gridinside=0.1;
  h=plot(d, opts); axis equal; title('test domain: exterior');end

% ext domain with holes
sq = segment.polyseglist(10, -1.7-0.5i + 0.5*[0 1i 1+1i 1]);  % CW square
sc = [segment(30, [1+1.5i, 0.5, 0, 2*pi/3]) segment(30, [1+1.5i, 0.5, 2*pi/3, 2*pi])]; % CCW circle via 2 segs
fprintf('exterior w/ 3 holes...\n'), disconnect(sa);
d = domain([], [], {sa sq sc}, {-1 ones(1,4) [-1 -1]});
if verb, figure; opts.gridinside=0.1;
  h=plot(d, opts); axis equal; title('test domain: exterior w/ holes'); end
