function u = unique(c)
% UNIQUE - remove duplicates from cell array of arbitrary objects
%
%  u = UNIQUE(c) where c is a cell array returns a cell row array u which has
%   all duplicates removed. The elements of c may be arbitrary class instances,
%   including user-defined classes.
%
%  Notes:
%   1) algorithm is currently O(N^2) since no known way to fast sort arbitrary
%   class instances. A fix would be to use sort() on doubles which uniquely
%   identify various class instances, as with graphics handles. 
%   2) this is needed since built-in unique only handles sort-able cell arrays,
%   ie string cell arrays

% alex barnett 8/13/08

if isempty(c), u = {}; return; end

u = {c{1}};
for i=2:numel(c)
  t = c{i}; tc = class(t);
  isnew = 1;
  for j=1:numel(u)
    if isa(u{j}, tc) & size(u{j})==size(t) & isequal(u{j}, t), isnew = 0; end
  end
  if isnew, u = {u{:} t}; end
end
