function circle(cx, cy, r, linetype);

N = 40;
x = zeros(1,N+1);
y = zeros(1,N+1);
for n=1:N+1
  x(n) = cx + r*cos(2*pi*n/N);
  y(n) = cy + r*sin(2*pi*n/N);
end
hold on; plot(x, y, linetype); hold off;
