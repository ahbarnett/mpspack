function addepwbasis(d, varargin)
% ADDEPWBASIS - create a real plane wave basis object in a domain
%
%  ADDEPWBASIS(d, N, opts) creates a real-valued plane wave basis
%   object within a domain object whose handle is d.
%   The rest of the argument list is discussed in RPWBASIS
%
% See also: RPWBASIS

% Copyright (C) 2008, 2009, Alex Barnett, Timo Betcke


d.bas  = {d.bas{:}, epwbasis(varargin{:})}; % append cell arr of basis handles

d.bas{end}.doms = d;                    % tell this basis it affects this domain
