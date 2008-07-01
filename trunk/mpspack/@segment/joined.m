% JOINED - # ways two segments have an endpoint in common (within eps)

function j = joined(s,t)
eps = 1e-15;          % corner locking distance
j = numel(find(abs(s.e(1)-t.e)<eps)) + numel(find(abs(s.e(2)-t.e)<eps));
