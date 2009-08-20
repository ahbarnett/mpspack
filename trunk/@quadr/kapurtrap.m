function [x,w,cs,ier]=kapurtrap(n,m)
%function [x,w,cs,ier]=kapurtrap(n,m)
%
% Construct the nodes and weights of n-point corrected trapezoidal 
% quadrature formula on the interval [0,1]. 
% Kapur-Rokhlin's correction formula.
%
% Periodic functions only!
%
% Input parameters:
%
% n - the number of quadrature points 
% m - the order of correction, valid values 2, 6, 10
%   
% Output parameters:
%
% n - integration nodes
% w - integration weights
% cs = end-point correction weights
%
% ier - the error code 
%         ier = 0    normal execution
%         ier = 8    m is not 2, 6, 10
%         ier = 16   n < m
%

% Copyright (C) 2008, 2009, Alex Barnett, Timo Betcke

ier=0;

T{1} = [1.825748064736159 -1.325748064736159];
T{2} = [4.967362978287758 -16.20501504859126 25.85153761832639 ...
        -22.22599466791883 9.930104998037539 -1.817995878141594];
T{3} = [7.832432020568779 -4.565161670374749e+1 1.452168846354677e+2 ...
        -2.901348302886379e+2 3.870862162579900e+2 -3.523821383570681e+2 ...
        2.172421547519342e+2 -8.707796087382991e+1 2.053584266072635e+1 ...
       -2.166984103403823];


x=(0:(n-1))'/(n-1);
w=ones(n,1);
%%w(1)=1/2; 
%%w(n)=1/2;
w(1)=0; % punctured trapezoidal formula 
w(n)=0; % punctured trapezoidal formula 

if( ~(m == 2 || m == 6 || m == 10 ) )
   ier = 8;
   w=w/(n-1);
   return
end

if( m > n )
   ier = 16;
   w=w/(n-1);
   return
end

if( m == 2 ) cs = T{1}; end
if( m == 6 ) cs = T{2}; end
if( m == 10 ) cs = T{3}; end


k=length(cs);
for j=1:k
    w(j+1)=w(j+1)+cs(j);
    w(n-j+1-1)=w(n-j+1-1)+cs(j);
end

w=w/(n-1);

