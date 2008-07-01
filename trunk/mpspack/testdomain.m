% test the domain constructor

clear all classes
p = [0 1 1+2i];
M = 10;
s(1) = segment(M,[p(1) p(2)]);
s(2) = segment(M,[p(2) p(3)]);
s(3) = segment(M,[p(1) p(3)]);
d = domain(s, [1 1 -1]);
figure; plot(d); axis equal; title('test domain: triangle');
[xx,yy] = meshgrid(-1:0.1:2,-1:0.1:2);
i = d.inside(xx+1i*yy);
hold on; plot(xx(i), yy(i), '.', 'markersize', 1);

if 0
  s(3) = segment(M,[p(1) p(3)]);
  d = domain(s,[1 1 -1]);
  d = domain(s,ones(1,length(s)));   % should cause corner error
end





