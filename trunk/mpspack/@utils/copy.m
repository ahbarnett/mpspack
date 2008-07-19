% COPY - Make a deep copy of a handle object
%
%  b = copy(a) supposedly makes b a deep copy of handle object a
%
%  author: Doug M. Schwarz, 6/16/08
%  http://www.mathworks.com/matlabcentral/newsreader/view_thread/171019#438411
% Also see
%  http://www.mathworks.com/matlabcentral/newsreader/view_thread/172147

function new = copy(this)

% Instantiate new object of the same class.
new = feval(class(this));

% Copy all non-hidden properties.
p = properties(this);
for i = 1:length(p)
  new.(p{i}) = this.(p{i});
end
