function [u xx yy ubwave] = showblochmode(sys, b, uc, s, vec, o, opts)
% temporary routine to by-hand plot Bloch modes from their singular vector vec
% Does empty, D/N metallic, dielectric. Will be superceded by @blochodeproblem
%
% [u xx yy ubwave] = showblochmode(sys, b, uc, s, vec, o, opts)
% 
% opts contains eval/plotting options:
% opts.bwcan - if true, cancels Bloch wave to leave periodic (non-PDE-soln) func
% opts.uctrim - if true, trims to UC (default true)
% opts.{gx,gy} - overrides x,y plot grids (either skewed or not)
% opts.skew - if true, interprets gx,gy as grids in UC coords (UC is [0,1]^2)
% opts.noplot - don't plot anything, just evaluate
% opts.Jfilter - if present, pass to QP field eval ('bo')
% opts.fixwrap - if true, attempt to correctly wrap field vals, only right (1,0)
%
% Example usage: see testblochbyhand.m, fig_blochmode.m
%
% Barnett 1/27/08, skewed 2/3/09, wrap field value from right hack 1/25/10

if nargin<7, opts = []; end
if ~isfield(opts, 'bwcan'), opts.bwcan = 0; end
if ~isfield(opts, 'uctrim'), opts.uctrim = 1; end, uctrim = opts.uctrim;
skew = isfield(opts,'skew');
if skew, N = 100; gx = ((1:N)-0.5)/N - 0.5; gy = gx;  % default UC-skewed plot
else, gx = -.75:0.01:.75; gy = -.5:0.01:.5; end       % default Cartesian plot
if isfield(opts, 'gx'), gx = opts.gx; end
if isfield(opts, 'gy'), gy = opts.gy; end
if ~isfield(opts, 'fixwrap'), opts.fixwrap = 0; end

[xx yy] = meshgrid(gx,gy);        % usual Cartesian grid
if isfield(opts,'skew')
  T = [real(uc.e1) real(uc.e2); imag(uc.e1) imag(uc.e2)]; % UC transf matrix
  xy = T*[xx(:)'; yy(:)'];
  xx = reshape(xy(1,:), size(xx)); yy = reshape(xy(2,:), size(xx)); clear xy;
end
p = pointset(xx(:)+1i*yy(:));     % points to evaluate at
%if uctrim, p = pointset(p.x(find(uc.inside(p.x)))); end % keep only inside UC
if sys~='e'
  ii = s.dom{1}.inside(p.x);          % boolean for pts `inside' exterior domain
  if opts.fixwrap, iir = s.dom{2}.inside(p.x+uc.e1); ii = ii & ~iir; %(1,0) only
  else iir = find(0*xx); end  % iir true if x is in incl neighbor to left
  if uctrim, ii = ii & uc.inside(p.x); end           % keep only inside UC
  po = pointset(p.x(ii));
else, ii = find(1+0*xx); po = p; end  % indices of all pts

N = numel(vec);

disp('eval QP field...')
bo = []; if isfield(opts,'Jfilter'), bo.Jfilter=opts.Jfilter; end
tic; B = uc.evalbases(po, bo); fprintf('\tdone in %g s\n', toc) % QPLP eval mat
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
  tic; u(ii) = u(ii) + ext_dlp.evalunitcellcopies(po, uc, o) * c_dlp;
  fprintf('\tdone in %g s\n', toc)
  o.dom = s.dom{2};                         % now evaluate in interior domain
  ii = ~ii; if uctrim, ii = ii & uc.inside(p.x); end  % keep only inside UC
  x = p.x; if opts.fixwrap, x = x+uc.e1*iir; end % shift points rightwards
  po = pointset(x(ii));                       % ~ii = Boolean for int obst pts
  disp('eval interior obst field...');      % note no copies summation here
  tic; u(ii) = u(ii) + int_slp.eval(po)*c_slp + int_dlp.eval(po)*c_dlp;
  fprintf('\tdone in %g s\n', toc)
  if opts.fixwrap, u(iir) = u(iir) * uc.a^(-1); end % undo phase on shifted u
end

u = reshape(u, size(xx));
ubwave = u .* exp(-1i*real(conj(uc.kbloch).*(xx+1i*yy))); % conjugate bloch wave
if ~isfield(opts, 'noplot')
  figure;
  if opts.bwcan
    uplot = real(ubwave);
    ts = sprintf('k=%g, a=%gpi, b=%gpi: Re[u] (Bloch wave cancel)', ...
                 uc.k, angle(uc.a)/pi, angle(uc.b)/pi);
  else
    uplot = real(u);
    ts = sprintf('k=%g, a=%gpi, b=%gpi: Re[u]', uc.k, angle(uc.a)/pi, ...
                 angle(uc.b)/pi);
  end
  if skew, surf(xx, yy, 0*xx, uplot); shading interp; view(2); grid off;
  else imagesc(gx, gy, uplot); set(gca, 'ydir','normal');
  end
  op.arrow=0; op.label=0; op.normals=0; uc.plot(op);    % switch
  if sys~='e', s.plot(1, op); end
  xlabel x; ylabel y; axis equal; caxis(0.02*[-1 1]);
end
