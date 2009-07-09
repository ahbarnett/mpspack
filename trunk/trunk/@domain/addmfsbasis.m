% ADDMFSBASIS - create an MFS (fundamental solutions) basis object in a domain
%
%  ADDMFSBASIS(d, Z, tau, N, k, opts) creates a fundamental solutions basis
%   object within a domain object whose handle is d, with other arguments
%   as in MFSBASIS.
%   A warning is given if any charge point is inside the domain.
%
%  ADDMFSBASIS(d, [], tau, N, k, opts) makes an intelligent choice of MFS charge
%   points based on the domain. If the domain has one segment, it is therefore
%   closed, and MFS points will be placed in or out as appropriate. If there is
%   more than onse segment, the option opts.segnum will choose which
%   disconnected piece is to be used to place the MFS points.
%   Currently only isolated-segment domains are
%   implemented, that is, each segment closes on itself. If tau is empty, a
%   default value will be chosen which places the points outside the domain.
%
% See also: MFSBASIS
%
% To do: include better methods to choose charge locations given things about
% the domain - this is hard in general! See our JCP 2008 paper.

function addmfsbasis(d, varargin)

if isempty(varargin{1})                     % Z was empty
  % automatically generate sensible charge curve for the domain...
  tau = varargin{2};
  if isempty(tau)
    tau = -d.pm * 0.1;                    % default tau, but outside domain
  end
  if numel(varargin)>4, opts = varargin{5}; end
  if numel(d.seg)==1                        % one segment, isolated, closed
    Z = @(t) d.seg(1).Z(t/2/pi);            % use the segment's parametric func
  elseif isfield(opts, 'segnum')            % >1 segment, but all isolated
    Z = @(t) d.seg(opts.segnum).Z(t/2/pi);  % use the segment's parametric func
  else
    fprintf('warning: addmfsbasis doesnt know a good way to choose MFS pts!\n')
    fprintf('Hint: if domain has only closed segments, choose w/ opts.segnum\n')
  end
  b = mfsbasis(Z, tau, varargin{3:end});    % pass in the above Z and tau
else
  b = mfsbasis(varargin{:});                % just pass arguments through
end

d.bas  = {d.bas{:}, b};                     % append cell arr of basis handles

if numel(varargin)>1
  d.k = varargin{2};                        % resets domain wavenumber
end

n_ins = numel(find(d.inside(b.y)));
if n_ins>0                                  % if any MFS pts inside domain
  fprintf('warning: %d MFS points are inside domain!\n', n_ins)
end
