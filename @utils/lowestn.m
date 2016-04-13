function y = lowestn(x,n)
% LOWESTN  return smallest n items in a list of numbers.
y = sort(x,'ascend');
y = y(1:n);
