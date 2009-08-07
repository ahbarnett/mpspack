function x = stackquadpts(s, pm)
% STACKQUADPTS - helper routine, ordered quad pts from signed connected seg list

% Copyright (C) 2008, 2009, Alex Barnett, Timo Betcke


x = [];
for j=1:numel(s)
  if pm(j)==1
    x = [x; s(j).x];
  else                            % reverse order since flipped seg
    x = [x; s(j).x(end:-1:1)];
  end
end
