function addmfsbasis(d, varargin)
% ADDMFSBASIS - Add an MFS basis to a domain
%
%  For a full list of possible arguments see MFSBASIS
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

