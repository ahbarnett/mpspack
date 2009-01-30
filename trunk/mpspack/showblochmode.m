function [u ubwave] = showblochmode(sys, b, uc, s, vec, o, opts)
% temporary routine to by-hand plot Bloch modes from their singular vector vec
% Does empty, D/N metallic, dielectric. Will be superceded by @blochodeproblem
%
% showblochmode(sys, b, uc, s, vec, o, opts)
% 
% opts.bwcan - if true, cancels Bloch wave to leave periodic (non-PDE-soln) func
% opts.uctrim - true trims to UC (default)
%
% Barnett 1/27/08

if nargin<7, opts = []; end
if ~isfield(opts, 'bwcan'), opts.bwcan = 0; end
if ~isfield(opts, 'uctrim'), opts.uctrim = 1; end, uctrim = opts.uctrim;

gx = -.75:0.01:.75; gy = -.5:0.01:.5;
[xx yy]=meshgrid(gx,gy); p=pointset(xx(:)+1i*yy(:)); 
if sys~='e'
  ii = s.dom{1}.inside(p.x);          % boolean for pts `inside' exterior domain
  if uctrim, ii = ii & uc.inside(p.x); end           % keep only inside UC
  po = pointset(p.x(ii));
else, ii = find(1+0*xx); po = p; end  % indices of all pts

N = numel(vec);

disp('eval QP field...')
tic; B = uc.evalbases(po); fprintf('\tdone in %g s\n', toc) % QPLP eval mat
Nq = size(B,2); No = N-Nq;            % # QP dofs, # obst dofs (0 if empty)
u = 0*xx; u(ii) = B*vec(No+1:end);                       % QP field

if sys=='s' | sys=='n' % ------- homog BC metal inclusion (n sets du/dn=0)
  o.dom = s.dom{1};             % evaluate in exterior domain
  disp('eval obst field...');
  A = b.evalunitcellcopies(po, uc, o); % obst eval mat
  fprintf('\tdone in %g s\n', toc)
  u(ii) = u(ii) + A*vec(1:No);                % u = QP field + obst field
 
elseif sys=='t'     % ---------- transmission dielectric incl (disconnected)
  ext_slp=b(1); ext_dlp=b(2); int_slp=b(3); int_dlp=b(4);  % basis handles
  c_dlp = vec(1:No/2);  c_slp = vec(No/2+1:No);    % coeff subvectors, NB order
  disp('eval ext SLP obst field...');
  u(ii) = u(ii) + ext_slp.evalunitcellcopies(po, uc, o) * c_slp;
  fprintf('\tdone in %g s\n', toc); disp('eval ext DLP obst field...');
  u(ii) = u(ii) + ext_dlp.evalunitcellcopies(po, uc, o) * c_dlp;
  fprintf('\tdone in %g s\n', toc)
  o.dom = s.dom{2};                         % now evaluate in interior domain
  ii = ~ii; if uctrim, ii = ii & uc.inside(p.x); end  % keep only inside UC
  po = pointset(p.x(ii));                  % ~ii = Boolean for int obst pts
  disp('eval interior obst field...');      % note no copies summation here
  u(ii) = u(ii) + int_slp.eval(po)*c_slp + int_dlp.eval(po)*c_dlp;
  fprintf('\tdone in %g s\n', toc)
end

u = reshape(u, size(xx));
ubwave = u .* exp(-1i*real(conj(uc.kbloch).*(xx+1i*yy))); % conjugate bloch wave
figure; if opts.bwcan, imagesc(gx, gy, real(ubwave));
title(sprintf('k=%g, a=%gpi, b=%gpi: Re[u] (Bloch wave cancel)', uc.k, angle(uc.a)/pi, angle(uc.b)/pi));
  else, imagesc(gx, gy, real(u)); title(sprintf('k=%g, a=%gpi, b=%gpi: Re[u]', uc.k, angle(uc.a)/pi, angle(uc.b)/pi));
end
set(gca, 'ydir','normal');
op.arrow=0; op.label=0; op.normals=0; uc.plot(op);    % switch
if sys~='e', s.plot(1, op); end
xlabel x; ylabel y; axis equal; caxis(0.02*[-1 1]);
