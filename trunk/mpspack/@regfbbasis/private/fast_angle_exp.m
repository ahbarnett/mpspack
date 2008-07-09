function [ang] = fast_angle_exp(n, z)
% function [ang] = fast_angle_exp(n, z)
%
% returns matrix of all rotations of unit vecs in directions given by complex
% #s z and multiples n. Ie if zh = z./abs(z) then return matrix of exp(1i*n*zh)
% n and z may be row or col vecs. ang has size n_z by n_n. n must be contiguous
% increasing and pass from negative to positive integers.
% Now I included 8-blocking of the loop to take advantage of matlab vectorized
% multiply better, about 2.5x faster.
%
% To do: case n(end)<-n(1) still untested! Could also do recursive blocking
%
% barnett 2/29/08, fixed NaN when z=0 4/3/08

M = numel(z); z = reshape(z, [M 1]);
N = numel(n); n = reshape(n, [1 N]);
noff = 1-n(1);                                  % n index offset
nun = 8;                                        % block size
ang = zeros(M, N);
ang(:,noff) = 1;
z(find(z==0)) = 1;                              % dummy vals replace case z=0
z = z ./ abs(z);                                % put on unit circle

if n(end)>=-n(1)                % ............... do pos n first, copy for neg
  nst = 1;                                      % n to start at
  if n(end)>nun
    for l=1:nun-1                               % fill the first block
      ang(:,noff+l) = ang(:,noff+l-1) .* z;
    end
    znun = z.^nun;
    for b=1:floor(n(end)/nun - 1)               % may execute 0 times
      ang(:,noff+nun*b+(0:nun-1)) = ang(:,noff+nun*(b-1)+(0:nun-1)) .* ...
          repmat(znun, [1 nun]);
    end
    nst = floor(n(end)/nun) * nun;
  end
  for l=nst:n(end)                              % finish remaining ones naively
    ang(:,noff+l) = ang(:,noff+l-1) .* z;
  end
  ang(:,1:noff-1) = conj(ang(:,noff-n(1):-1:noff+1)); % copy only needed neg n

else                             % .............. do neg n first, copy for pos
  z = conj(z);                                  % since going down in n
  nst = -1;                                     % n to start at
  if -n(1)>nun
    for l=-1:-1:-nun+1                          % fill the first block
      ang(:,noff+l) = ang(:,noff+l+1) .* z;
    end
    znun = z.^nun;
    for b=1:floor(-n(1)/nun - 1)                % may execute 0 times
      ang(:,noff-nun*b-(0:nun-1)) = ang(:,noff-nun*(b-1)-(0:nun-1)) .* ...
          repmat(znun, [1 nun]);
    end
    nst = -floor(-n(1)/nun) * nun;
  end
  for l=nst:-1:n(1)
    ang(:,noff+l) = ang(:,noff+l+1) .* z;
  end
  ang(:,noff+1:end) = conj(ang(:,noff-1:-1:noff-n(end)));
end
% note the 2nd half code is simply the first half replacing n(end) by -n(1)
% and doing loops with a -1 step. Could probably be more elegantly written
