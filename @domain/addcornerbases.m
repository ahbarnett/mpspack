function addcornerbases(d, N, opts)
% ADDCORNERBASES - add irreg Fourier-Bessel basis to each corner of a domain
%
%  ADDCORNERBASES(d, N) adds one irregular fractional-order wedge FB basis
%    set at each corner of domain handle d, of total degree N (total # degrees
%    of freedom is 2*N). Degree of each corner is proportional to angle.
%
%  ADDCORNERBASES(d, N, opts) is as above, but allows options to be chosen:
%    opts.cornermultipliers = nmultipliers (basis size proportions) for each
%                       corner. If zero for a corner, no basis is added there.
%                       (default is all equal to 1 for valid corners)
%    opts.type = 's', 'c', or 'cs', chooses angular sin/cos type as in nufbbasis
%                       (default is 'cs')
%    Other opts fields are passed to nufbbasis.
%
% See also: NUFBBASIS, ADDNUFBBASIS
%
% Copyright (C) 2008, 2009, Timo Betcke, Alex Barnett


  if nargin<3, opts = []; end
  if ~isfield(opts, 'cornermultipliers')   % default: 1 for each valid corner...
    opts.cornermultipliers = ones(size(d.cloc)) .* ~isnan(d.cloc);
  end
  %  ...might be better with d.cang/pi instead of 1 ?
  
  for j=1:numel(d.cloc)
    bra = -d.cangoff(j) - d.cang(j)/2;  % choose branch cut facing away
    if opts.cornermultipliers(j)>0
      opts.nmultiplier = opts.cornermultipliers(j);
      d.addnufbbasis(d.cloc(j), pi/d.cang(j), d.cangoff(j), bra, ...
                     ceil(N*d.cang(j)/pi), opts);
    end
  end

