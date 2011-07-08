function s = smoothfourier(M, aj, bj)
% SMOOTHFOURIER - general Fourier series radial function closed segment
%
%  s = SMOOTHFOURIER(M, aj, bj) generates a smooth closed radial function
%   segment with M discretization pts, Fourier cos amplitudes aj and Fourier
%   sin amplitudes bj. The zero-freq term is 1, and aj and bj are the coeffs
%   of indices 1...n. aj and bj may be row or col vectors, but numel(bj)
%   must be at least numel(aj)

% Copyright (C) 2008 - 2011, Alex Barnett, Timo Betcke

n = numel(aj);
if numel(bj)~=n, error('bj must have n=numel(aj) elements!'); end
aj = aj(:).'; bj = bj(:).';  % make row vecs
j = (1:n).';                 % col vec of indices
if n==0, s = segment.radialfunc(M, {@(t) 1+0*t, @(t) 0*t, @(t) 0*t});
else % now use outer products inside the trig funcs...
  s = segment.radialfunc(M, {@(t) 1 + reshape(aj*cos(j*t(:).') + bj*sin(j*t(:).'),size(t)), ...
                      @(t) reshape(-(j'.*aj)*sin(j*t(:).') +(j'.*bj)*cos(j*t(:).'),size(t)), ...
                      @(t) reshape(-(j'.^2.*aj)*cos(j*t(:).') - (j'.^2.*bj)*sin(j*t(:).') ,size(t))});
end
