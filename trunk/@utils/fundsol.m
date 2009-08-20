function [F0 F1 F2] = fundsol(r, k,orders, opts)
% FUNDSOL - Compute fundamental solutions
%
%  [F0 F1 F2] = FUNDSOL(r, k, orders, fast) computes given a wavenumber k
%  a matrix of distances r the matrix of fundamental solutions F0, its
%  first r-derivative F1 and its second r-derivative F2.
%
%  If k>0 then F0=(1i/4)*besselh(0,k*r).
%  If k=0 then F0=-(1/2/pi)*log(r).
%
%  opts is a structure that can contain various options. Currently
%  supported is o.fast. If o.fast>0 use fast Hankel routines by Rokhlin and
%  Greengard.
%
%
%  The string orders can take the following values:
%
%    orders = '0'   - Evaluate only F0
%    orders = '1'   - Evaluate only F1
%    orders = '01'  - Evaluate F0 and F1
%    orders = '012' - Evaluate F0, F1 and F2
%    orders = '12'  - Evaluate F1 and F2
%

fast=opts.fast;
F0=[]; F1=[]; F2=[];

if abs(k)>0,
    % Helmholtz case
    if fast==2
        [F0 F1]=utils.greengardrokhlinhank106(k*r);
    elseif fast==1
        [F0 F1]=utils.greengardrokhlinhank103(k*r);
    else
        if strcmp(orders,'0'), 
            F0=besselh(0,k*r);
        elseif strcmp(orders,'1'),
            F1=besselh(1,k*r);
        else
            F0=besselh(0,k*r); F1=besselh(1,k*r);
        end
    end
    F0=(1i/4)*F0; F1=-k*(1i/4)*F1;
    if orders(end)=='2', F2=-k^2*F0-F1./r; end
    if orders(1)~='0', F0=[]; end % Delete F0 again since it is not asked
                                  % for (only important for order='12')
else
    % Laplace case
    switch orders,
        case '0', F0=-(1/2/pi)*log(r);
        case '1', F1=-(1/2/pi)./r;
        case '01', F0=-(1/2/pi)*log(r); F1=-(1/2/pi)./r;
        case '012', F0=-(1/2/pi)*log(r); F1=-(1/2/pi)./r; F2=-F1./r;
        case '12', F1=-(1/2/pi)./r; F2=-F1./r;
    end
end

            

        
            
