function bas = addinoutlayerpots(seg, varargin)
% ADDINOUTLAYERPOTS - add Rokhlin-style layer-potentials both sides of segment
%
%  ADDINOUTLAYERPOTS(segs, a) adds a mixture of SLP + DLP (coeffs given
%   by 1-by-2 array a) with densities lying on the segment handles segs, to
%   affect the domains lying *both* side of each segment. If any segment in the
%   list segs is not connected on both sides, an error results.
%
% See also LAYERPOT, SEGMENT, DOMAIN.ADDLAYERPOT

% Copyright (C) 2010, Alex Barnett, Timo Betcke

bas = {};
for s=seg
  d1 = s.dom{1}; d2 = s.dom{2};   % the connected domains
  if isempty(d1) | isempty(d2)
    error('each segment must be connected to domains on both sides!'); end
  b = layerpot(s, varargin{:});
  b.doms = [d1 d2];                     % array, not cell array (bkwd-compat)
  bas  = {bas{:}, b};                   % append cell arr of basis handles
  d1.bas = {d1.bas{:}, b};
  d2.bas = {d2.bas{:}, b};
end
