function i = isin(b, c)
% ISIN - returns true if first argument is in cell or array given by second arg
%
%   Note: works on arbitrary cell elements. Takes O(N) time. Now allows c to
%   be a non-cell array.
i = 0;
if isempty(b), return; end
if ~iscell(c), c = num2cell(c); end       % insure it's a cell object
for j=1:numel(c)
  t = c{j};
  if isequal(t, b), i=1; return; end
end
