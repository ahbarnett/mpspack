function [A,B,C] = para_fit(e, f)
% function [A,B,C] = para_fit(e, f)
%
% given e array of 3 ordinates, f array of corresp func values
% returns parameters in f = A + C(x-B)^2

x = e(2)-e(1); w = e(3)-e(2);
y = f(1)-f(2); z = f(3)-f(2);
C = (x*z + w*y)/(x*w*(e(3)-e(1))); % solve para coeffs: f = A + C(x-B)^2
B = (e(3)+e(2) - z/(w*C))/2;  
A = f(2) - C*(e(2)-B)^2;
