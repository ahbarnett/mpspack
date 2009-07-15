function i = isin(b, c)
% ISIN - returns true if first argument is in cell array given by second arg
%
%   Note: works on arbitrary cell elements. Takes O(N) time.
i = 0;
if isempty(b), return; end
if numel(c)==1 & ~iscell(c), c = {c}; end       % insure it's a cell object
for j=1:numel(c)
  t = c{j};
  if isequal(t, b), i=1; return; end
end
