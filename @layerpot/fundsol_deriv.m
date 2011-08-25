function [B radderivs] = fundsol_deriv(r, cosphi, k, radderivs)
% function [B {, radderivs}] = fundsol_deriv(r, cosphi, k {, radderivs})
%
% Returns matrix of (source)normal-deriv of fundamental soln
% given dist matrix r, angle matrix cosphi, and wavenumber k.
% Uses symmetry (since zero argument is fast) if r appears symm (based on diag).
% If requested, returns radderivs for later reuse as optional input argument
% if cosphi changes but nothing else does
% (the user is assumed to give correct radderivs of appropriate size, k, etc).
% Note: there is an overall minus sign due to it being the src n-deriv.
%

% Copyright (C) 2008, 2009, Alex Barnett and Timo Betcke

wantrad = nargout>1;

if nargin>3                                              % reuse radderivs
  B = radderivs .* cosphi;                               % note covers all k
  
else                                                     % compute from scratch
  if k==0                           % ........ laplace 
    if wantrad
      radderivs = (1/2/pi) ./ r;
      B = radderivs .* cosphi;
    else
      B = (1/2/pi) * cosphi ./ r;
    end
  else                              % ......... helmholtz
    % if self-interactions (square & diag=999), assume symm, do upper tri only
    if size(r,1)==size(r,2) & norm(diag(r)-999)<1e-14
      %disp(sprintf('self, diag(r)=%g', r(1,1)));
      B = triu(besselh(1, 1, k*triu(r,1)),1);
      B = B.' + B;
      B(diagind(B)) = besselh(1, 1, k*diag(r));  % always dummy
    else  % do the usual thing which works for distant nonsymm interactions...
      B = besselh(1, 1, k*r);
    end
    % currently B contains radderivs without the ik/4 prefactor
    if wantrad
      radderivs = (1i*k/4) * B;
      B = radderivs .* cosphi;
    else
      B = (1i*k/4) * B .* cosphi;
    end
  end
end
