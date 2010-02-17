function [M data] = fillblochmodematrix(sys, b, uc, s, o, bo, opts)
% temporary routine to do by-hand filling of single-segment systems
% under test in testblockbyhand. (no early expts such as 3xUC, Dirichlet wrap)
% Does empty, D/N metallic, dielectric. Will be superceded by @blochmodeproblem
%
% sys = system, b = array of basis objects, uc = unit cell, s = closed segment
% o = opts for A,C,Q; bo = opts for B; opts = general opts for this routine
% Eg opts.data contains A,B,C data as data.{dA,dB,dC}
% bo.close : close evaluation (recommended only for B).
%
% Notes:
%  1) sys='n' Neumann fails with BWLP since no cot() hypersingular term yet
%   coded in @layerpot/T.m
% Barnett 1/27/08

if nargin<7, opts = []; end
if ~isfield(opts, 'verb'), opts.verb = 0; end, v = opts.verb;   % verbosity
if ~isfield(opts, 'wrap'), opts.wrap = 0; end, wr = opts.wrap;  % wrap flag
isdata = isfield(opts, 'data'); if isdata, data = opts.data; end

tic; Q = uc.evalbasesdiscrep(o);   % poly data storage kept in uc automatically
if v, fprintf('sys=%s\tfill Q: %g s\n', sys, toc); end

if sys=='e'             % ---------- empty unit-cell case
  M = Q; data = [];             % dummy data output

elseif sys=='s' | sys=='n' % ------- homog BC metal inclusion (n sets du/dn=0)
  if isdata, bo.data = data.dB; end
  tic; if sys=='s', [B data.dB] = uc.evalbaseswithdata(s, bo); % heeds Jfilter
       else [dummy B data.dB] = uc.evalbaseswithdata(s, bo); % use Bn
  end
  if v, fprintf('\tfill B: %g s\n', toc); end
  o.dom = s.dom{1};             % evaluate in exterior domain, for filling A...
  if isdata, o.data = data.dA; end
  tic; if sys=='s', [A data.dA] = b.evalunitcellcopies(s, uc, o); % sq src blk
      else [dummy A data.dA] = b.evalunitcellcopies(s, uc, o); % use An
  end
  if v, fprintf('\tfill A: %g s\n', toc); end
  if isdata, o.data = data.dC; end
  tic; [C data.dC] = b.evalunitcelldiscrep(uc, o); % heeds ucbuf
  if v, fprintf('\tfill C: %g s\n', toc); end
  M = [2*A 2*B; C Q];

  elseif sys=='t'     % ---------- transmission dielectric incl (disconnected)
  ext_slp=b(1); ext_dlp=b(2); int_slp=b(3); int_dlp=b(4);  % basis handles
  if isdata, bo.data = data.dB; end
  tic; [B Bn data.dB] = uc.evalbaseswithdata(s, bo);
  B = [B; Bn];
  if v, fprintf('\tfill B: %g s\n', toc); end
  tic; o.dom = s.dom{1};         % filling A: first eval in exterior domain...
  if isdata, o.data = data.dADo; end
  [ADo ATo data.dADo] = ext_dlp.evalunitcellcopies(s, uc, o);
  if isdata, o.data = data.dASo; end
  [ASo ADTo data.dASo] = ext_slp.evalunitcellcopies(s, uc, o);
  oi = o; oi.dom = s.dom{2}; oi.nei=0;  % now use int dom, no copies (nei=0)...
  if isdata, oi.data = data.dADi; end   % (uccopies used since does a-poly too)
  [ADi ATi data.dADi] = int_dlp.evalunitcellcopies(s, uc, oi);  % NB oi opts
  if isdata, oi.data = data.dASi; end
  [ASi ADTi data.dASi] = int_slp.evalunitcellcopies(s, uc, oi);
  A = [ADo-ADi, ASo-ASi; ATo-ATi, ADTo-ADTi]; % diel mism, 2x2 block
  if v, fprintf('\tfill A: %g s\n', toc); end
  tic; if isdata, o.data = data.dCD; end
  [CD data.dCD] = ext_dlp.evalunitcelldiscrep(uc, o);
  if isdata, o.data = data.dCS; end
  [CS data.dCS] = ext_slp.evalunitcelldiscrep(uc, o);
  C = [CD CS];                   % coeff subvector order: obst DLP before SLP
  if v, fprintf('\tfill C: %g s\n', toc); end
  M = [2*A 2*B; C Q];
end

