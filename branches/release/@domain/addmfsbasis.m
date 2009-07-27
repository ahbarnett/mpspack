function addmfsbasis(d, varargin)
% ADDMFSBASIS - Add an MFS (fundamental solutions) basis to a domain
%
% d.ADDMFSBASIS(args) where d is a domain, adds an MFS basis to the domain.
%   The argument list args are exactly as in MFSBASIS.
%
% See also: MFSBASIS

% Copyright (C) 2008, 2009, Timo Betcke, Alex Barnett

b=mfsbasis(varargin{:});
d.bas={d.bas{:},b};                    % append basis to domain's list
d.bas{end}.doms = d;                   % tell this basis it affects this domain
n_ins = numel(find(d.inside(b.y)));
if n_ins>0                                  % if any MFS pts inside domain
  fprintf('warning: %d MFS points are inside domain!\n', n_ins)
end

