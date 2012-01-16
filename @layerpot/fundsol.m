function A = fundsol(r, k)
% function A = fundsol(r, k)
%
% Returns matrix of fundamental solns A given dist matrix r, and wavenumber k.
% Uses Matlab built-in hankels, and
% uses symmetry (since zero argument is fast) if r appears symm (based on
% flagging all diag vals with an impossible flag value).
%
% See also LAYERPOT

% Copyright (C) 2008 - 2012, Alex Barnett and Timo Betcke

symmflagval = -999;   % all diag vals of this signifies symmetric - a hack

%fprintf('min r = %g\n', min(r(:)))
if k==0
  A = -(1/2/pi) * log(r);                                % laplace 
else
  % if self-interactions (square & dummy diag), assume symm, do upper tri only
  if size(r,1)==size(r,2) & norm(diag(r)-symmflagval)<1e-14  % hack!
    %disp(sprintf('fundsol symm: diag(r)=%g, N=%d', r(1,1), size(r,1)));
    A = (1i/4) * triu(besselh(0, 1, k*triu(r,1)),1);     % helmholtz
    A = A.' + A;
    A(diagind(A)) = (1i/4) * besselh(0, 1, k*diag(r));
  else  % do the usual thing which works for distant nonsymm interactions...
        %disp(sprintf('fundsol unsymm: r(1,1)=%g, M=%d, N=%d', r(1,1), size(r,1),...
        %size(r,2)));
    A = (1i/4) * besselh(0, 1, k*r);                     % helmholtz
  end
end
