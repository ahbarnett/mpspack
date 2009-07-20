function addcornerbases(d, N, opts)
% ADDCORNERBASES - add irreg Fourier-Bessel basis to each corner of a domain
%
%  ADDCORNERBASES(d, N) adds one irregular fractional-order wedge FB basis
%    set at each corner of domain handle d, of total degree N (total # degrees
%    of freedom is 2*N). Degree of each corner is proportional to angle.
%
%  ADDCORNERBASES(d, N, opts) is as above, but allows options to be chosen:
%    opts.cornerflags = list of numbers, either false or true, specifying which
%                       corners are to be added (default is all true)
%    opts.type = 's', 'c', or 'cs', chooses angular sin/cos type as in nufbbasis
%                       (default is 'cs')
%    Other opts fields are passed to nufbbasis.
%
% See also: NUFBBASIS, ADDNUFBBASIS
%
% Copyright (C) 2008, 2009, Timo Betcke, Alex Barnett


  if nargin<3, opts = []; end
  if ~isfield(opts, 'cornerflags'), opts.cornerflags = ones(size(d.cloc)); end
  
  for j=1:numel(d.cloc)
    bra = -d.cangoff(j) - d.cang(j)/2;  % choose branch cut facing away
    if opts.cornerflags(j)
      d.addnufbbasis(d.cloc(j), pi/d.cang(j), d.cangoff(j), bra, ...
                     ceil(N*d.cang(j)/pi), opts);
      d.bas{end}.nmultiplier = d.cang(j)/pi;    % hack idea now - TIMO FIX
    end
  end

