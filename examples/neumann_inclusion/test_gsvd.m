% understanding GSVD and it's two forms of X output. Matlab vs G-vL, etc.
% Barnett 12/7/15
clear
N = 1e2;
A = rand(N);
B = rand(N);
[UU VV X C S] = gsvd(A,B);
t = sqrt(diag(C'*C)./diag(S'*S));
X = inv(X');   % crucial, since Matlab's X defined differently!
l = find(t==min(t)); t = t(l); x = X(:,l);  % l=index of min gsingval

t
norm(A*x)/norm(B*x)
% are equal
