function [G An] = intnormmatrix(s, k, A, Dx, Dy)
% G = intnormmatrix(d), s is a segment, A is basis value matrix, Dx, Dy the
%  basis deriv matrices (rectangular).
% G = interior norm matrix quadratic form, PSD
% [G An] = ... returns normal-deriv matrix
%
% For use by tensionsq & tension - move into github when done
% Barnett 12/7/15

N = size(A,2);                      % # basis funcs (or rotated basis)
xn = real(conj(s.x) .* s.nx);        % x.n weight factor on bdry
E = k^2;

An = repmat(real(s.nx), [1 N]).*Dx + repmat(imag(s.nx), [1 N]).*Dy;
W = repmat(s.w.'.*xn, [1 N]);              % x.n as elementwise mult to right
G = A' * (W .* A) * E;
G = G - Dx' * (W .* Dx) - Dy' * (W .* Dy);
W = repmat(s.w.', [1 N]);                  % 1 as elementwise mult to right
Ddil = repmat(real(s.x), [1 N]).*Dx + repmat(imag(s.x), [1 N]).*Dy;
B = Ddil'*(W.*An);
G = G + B + B';                         % conjugate saves another mat mult
G = (G+G')/(4*E);                       % rescale and make symm
