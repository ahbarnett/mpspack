% test constructors, how slow, eg for segments
% barnett 8/11/08

N = 10; n = 300;
tic;
for i=1:n
  s(i) = segment(N, [0 1+1i]);
end
fprintf('segment create: %.2g ms\n', toc*1000/n);
tic;
for i=1:n  % note may exceed recursion depth on s.Z etc func handles
  s(i).translate(1-2i);
end
fprintf('segment translate: %.2g ms\n', toc*1000/n);
tic;
%profile clear; profile on;
for i=1:n
  t = s(i).translate(1-2i);
end
%profile off; profile viewer
fprintf('segment translate making new: %.2g ms\n', toc*1000/n);
tic;
for i=1:n
  lp = layerpot(s(i), 's', 10);
end
fprintf('layerpot create: %.2g ms\n', toc*1000/n);
c = cell(1, n); for i=1:n, c{i} = s(i); end
tic; utils.isin(s(floor(n/2)), c);
fprintf('utils.isin (isequal): %.2g ms\n', toc*1000/n);
