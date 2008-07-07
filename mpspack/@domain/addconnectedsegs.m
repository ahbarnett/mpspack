% ADDCONNECTEDSEGMENTS - append a domain's params, corners given conn seg list 
%
%  This is a helper routine for domain constructor. Doc to be written

function d = addconnectedsegs(d, s, pm, o)

if nargin<4, o = []; end
if ~isfield(o, 'hole'), o.hole = 0; end    % default is an outer bdry (+ve area)

if isempty(s), return; end                 % nothing to do if no segments added
s = reshape(s, [1 numel(s)]); pm = reshape(pm, [1 numel(pm)]); 
d.seg = [d.seg s]; d.pm = [d.pm pm];       % append segments and signs as is
if ~o.hole
  d.spiece = zeros(size(s));               % outer bdry always piece # 0 (pri)
else
  if isempty(d.spiece), nextp = 1; else nextp = d.spiece(end)+1; end
  d.spiece = [d.spiece, nextp * ones(size(s))];  % next piece label
end
d.perim = d.perim + sum([s.w]);            % sum of all quadr weights
allx = domain.stackquadpts(s, pm);         % get quad pts in s in correct order
a = polyarea(real(allx), imag(allx));      % area: use all quadr pts as polygon
if o.hole
  d.area = d.area - a;
else
  d.area = d.area + a;
end

% create corners...
% indices (1 or 2) of which end starts the segment ahead of each corner...
i = (1-pm)/2+1;
% now make versions of the above 3 lists, for the previous segment...
prevs = circshift(s, [0 1]); prevpm = circshift(pm, [0 1]);
previ = (1+prevpm)/2+1; % ind (1 or 2) for segment behind each corner
eps = 1e-15;                     % max allowed corner position error
for j=1:length(s)                          % loop over possible corners in list
  d.cloc = [d.cloc s(j).eloc(i(j))];       % append to corner loc list
  if abs(d.cloc(end) - prevs(j).eloc(previ(j))) > eps
    fprintf('domain warning: corner %d not conn (normalscheck invalid)!\n', j);
    d.cloc(end) = NaN;                     % invalidate this corner
  end
  d.cangoff = [d.cangoff, pm(j)*s(j).eang(i(j))]; % append, sign flipped by pm
  cang = -prevpm(j)*prevs(j).eang(previ(j))/d.cangoff(end); % on unit circle
  % now need to take the log with branch cut along +ve real axis...
  d.cang = [d.cang, pi + imag(log(-cang))]; % this is angle in (0,2pi)
                                 % note nothing changes here if hole domain
end
