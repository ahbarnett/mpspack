% DATAWRAPR - trial wrapping-over-R-wall of mat data, eg a cell item in dB
%
% used by testblochwrap.m in the 1xUC scheme. Moves rows of each matrix in d.B
% to one higher alpha-power. (Tries not to create too many new
% entries in copylist? - no, just doubles everything... fine for now!)
%
% d = uc.datawrapR(d, i)   d is data struc, i is list of indices to wrap
function d = datawrapR(uc, d, i)
  c = d.copylist;
  c.apow = [c.apow c.apow+1];       % increase alpha power by 1
  c.bpow = [c.bpow c.bpow]; c.t = [c.t c.t]; c.remph = [c.remph c.remph];
  n = size(d.B,3);       % original size of copylist
  d.B(:,:,n+(1:n)) = 0;         % new empty matrix blocks
  d.B(i,:,n+(1:n)) = d.B(i,:,1:n);    % rows to move
  d.B(i,:,1:n) = 0;                   % kill the rows which moved out
  d.copylist = c;
end

% More ambitious sorting code, not finished, would create less data...
%  [c.apow j] = sort(c.apow, 'descend');   % work from largest to smallest a pows
%  c.bpow c.bpow(j); c.t = c.t(j); c.remph = c.remph(j);  % shuffle
%  no = size(d.B,3);       % original size of copylist
%  n = no;                 % current size
%  for k=1:n
%    ap = c.apow(k)+1;        % new apow
%    if ~isempty(find(c.apow+1==))  % sth already at that (a,b)-pow, so add to it
%      
%    else
%      
%    end
%  end
%  d.copylist = c;
%end