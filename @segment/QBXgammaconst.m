function gam = QBXgammaconst(s, t, R, opts)
%
% gam = QBXgammaconst(seg, t, R) returns the gamma constant for the
%  segment seg, at the complex [0,1]-scaled parameter location
%  t, for radius R. (R is just a multiplicative factor.) It does this by
%  crudely sampling distances to a range of complexified segments.
%
% Barnett 9/25/12

n = numel(s.t) * 10; % napproxv
nd = 20; % how many distances to sample
al0 = imag(t);
maxaoverd = 0;            % max alpha-dist over actual dist
for i=1:nd, al = (i-1)/nd*al0;  % imag dist, for [0,1] scaled seg
  dist = min(abs(s.Z(t) - s.Z((1:n)/n + 1i*al)));
  maxaoverd = max(maxaoverd, (al0-al)/dist);
end
gam = maxaoverd * R / al0;
