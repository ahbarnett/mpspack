function [a b] = dielectriccoeffs(pol, np, nm) % ...diel matching coeffs
% DIELECTRICCOEFFS - give a and b coeff pairs (1-by-2) from refractive indices
%
%  [a b] = dielectriccoeffs(pol, np, nm). np and nm are indices on + and -
%   sides respectively
%
% Copyright (C) 2008, 2009, Timo Betcke, Alex Barnett


epsp = np^2; epsm = nm^2;               % epsilon (permittivity)
if pol=='tm' | pol=='TM'                % u represents Ez
  a = [1 -1]; b = [1 -1];               % note opposing signs, for continuity!
elseif pol=='te' | pol=='TE'            % u represents Hz
  a = [1 -1]; b = [1/epsp -1/epsm];     % 1/n^2 u_n is continuous (Wiersig '02)
else
  error('polarization must be TM or TE');
end
