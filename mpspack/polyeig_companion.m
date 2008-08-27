function [X E] = polyeig_companion(varargin)
% hack simple first companion matrix form for PEP (Mehrmann's 2004 article)
% barnett 8/20/08
  k = numel(varargin)-1;
  n = size(varargin{1},1);
  A = -eye(n*k); B = diag(-ones(1,n*(k-1)), n);
  A(1:n,1:n) = -varargin{k+1};
  for i=0:k-1
    B((1:n)+(k-1-i)*n,1:n) = varargin{i+1};
  end
  if nargout==1
    X = eig(B,A);
  else
    [X,E] = eig(B,A);
    X = X(1:n,:);
  end
  