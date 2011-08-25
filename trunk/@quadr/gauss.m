  function [x,w] = gauss(N)
% GAUSS  nodes x (Legendre points) and weights w
%        for Gauss quadrature

% Copyright (C) 2008, 2009, Alex Barnett, Timo Betcke
  persistent xstore wstore Nstore
  if N==Nstore,
      x=xstore; w=wstore;
      return
  end
  if N>300
    fprintf('warning: finding >300 gauss quadr pts slow O(M^3)!\n');
  end
  beta = .5./sqrt(1-(2*(1:N-1)).^(-2));
  T = diag(beta,1) + diag(beta,-1);
  [V,D] = eig(T);
  x = diag(D); [x,i] = sort(x);
  w = 2*V(1,i).^2;
  Nstore=N; wstore=w; xstore=x;
