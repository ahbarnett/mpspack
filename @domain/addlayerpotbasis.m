function bas = addlayerpotbasis(d, seg, varargin)
% ADDLAYERPOTBASIS - create a layer-potential basis set object in a domain
%
%  ADDLAYERPOTBASIS(d, segs, a) adds a mixture of SLP + DLP (coeffs given
%   in 1-by-2 array a) with densities on the segment handles segs, to domain d.
%   Note that segs may be segments not on the boundary of d.
%   If segs is empty, the boundary of d will be assumed. The mixture is fixed
%   across segments; if you wish to choose different mixtures, use multiple
%   calls to addlayerpotbasis.
%
%  ADDLAYERPOTBASIS(d, segs, a, opts) passes options, see opts in LAYERPOT.
%
%  bas = ADDLAYERPOTBASIS(...) returns cell array of handle(s) of created
%   basis object(s).
%
% See also LAYERPOT
%
%  Notes/issues:
%  To do: change to add images, or multiple segs as one basis object?

if isempty(seg), seg = d.seg; end
bas = {};
for s=seg
  b = layerpot(s, varargin{:});
  bas  = {bas{:}, b};                     % append cell arr of basis handles
end
d.bas = {d.bas{:}, bas{:}};               % append to domain's referred bases
