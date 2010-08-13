% EVP - define a Laplace eigenvalue problem for a single interior domain
%
% p = evp(dom) creates an EVP object from the domain dom, comprising the
%   domain along with its segments the boundary conditions they carry, and
%   the basis sets associated with the domain.
%
%   Preliminary version: only one domain is supported, with Dirichlet BCs,
%   and the only solution method is 
%   root-search on the Fredholm determinant of Id minus double-layer operator.

% Copyright (C) 2010 Alex Barnett, Timo Betcke

classdef evp < problem & handle
  properties
    kj                        % eigenwavenumbers = sqrt(eigenvalues) (col vec)
    coj                       % basis coeffs (columns of array) of modes
    ej                        % estimated errors of kj (col vec)
    ndj                       % normal derivatives (cols of array) of modes
  end

  methods % ----------------------------------------- particular to EVP
    function pr = evp(doms) % .........................constructor
      if nargin==0, return; end    % needs empty constructor
      if numel(doms)>1, error('only single domain currently supported!'); end
      pr.doms = doms; % Following code duplicates bvp, should be farmed out:
      pr.segs = [];             % now build a list of segments, crude O(N^2)
      for d=doms
        for s=d.seg
          if isempty(find(pr.segs==s))  % only add the segment if it's new
            pr.segs = [pr.segs s];
          end
        end
      end
    end
  
    function [kj coj ndj ej] = solvespectrum(p, kwin, meth, o)
    % SOLVESPECTRUM - compute spectrum, and possibly modes, in wavenumber window
    %
    % [kj] = SOLVESPECTRUM(p, [klo khi]) uses a default method to find all
    %   eigenwavenumbers kj lying in the wavenumber interval [klo, khi], for the
    %   eigenvalue problem object p.
    %   Laplacian eigenvalues are simply the squares of the kj.
    %   The problem also stores the kj list as p.kj
    %
    % [kj coj ndj ej] = ... also returns basis coeffcients, normal derivatives,
    %   and error estimates, if possible.
    %
    % [...] = SOLVESPECTRUM(p, [klo khi], meth, opts) allows one to choose the
    %   method,
    %     meth = 'fd' : Fredholm det (layer-potential bases only)
    %                   Error estimates are dist to real axis of Boyd roots
    %   and choose options,
    %     opts.verb = 0,1,... : verbosity (0=silent, 1=diagnostic text, ...)
    %
    % Currently the only solution method is root-search on the Fredholm
    %  determinant of Id minus double-layer operator, for a single domain, with
    %  Dirichlet BCs everywhere on the boundary.
    %
    % Notes / to do: In development!
    % * keep more data of all the det evals, pass out?
      if nargin<3, meth='fd'; end   % default method
      if nargin<4, o = []; end
      if ~isfield(o, 'tol'), o.tol = 1e-8; end   % requested tolerance
      if ~isfield(o, 'verb'), o.verb = 1; end; v = o.verb; % default verbosity
      if numel(p.doms)~=1, error('problem must have exactly one domain!'); end
      d = p.doms;                            % the domain
      wantmodes = nargout>1;
      klo = kwin(1); khi = kwin(2);
      p.fillquadwei; M = numel(p.sqrtwei);   % total # boundary pts
      if v, fprintf('mean M ppw @ khi = %.3g\n', 2*pi*M/khi/d.perim); end
      
      if strcmp(meth,'fd') % ............... Fredholm det method
        k = klo; kj = []; ej = []; N = M;  % square system
        kscale = 70/d.diam;                % est k at which level density ~ 20
        io.Ftol = 1e-10; io.tol = o.tol;   % allows leeway in collecting roots
        io.disp = v;                       % opts for interval search
        while k<khi
          ktop = min(k + kscale*min(1/k,0.2), khi); % top of current k window
          % option to reset N here, requadrature all segs (in proportion?)
          % ...
          %f = @(k) det(eye(N) - 2*layerpot.D(k, s)); % one closed segment
          f = @(k) det(p.fillfredholmop(k)); % function to rootfind on
          if v, fprintf('rootfind in [%.16g, %.16g], N=%d:\n', k, ktop, N); end
          [kl el y u] = utils.intervalrootsboyd(f, [k ktop], io);
          [kl i] = sort(kl); el = el(i);
          if v, fprintf('found %d eigenvalues w/ %d fredholm det evals (ratio %.3g)\n', numel(kl), numel(y), numel(y)/numel(kl)); end
          kj = [kj; kl]; ej = [ej; el];
          k = ktop;                        % increment k and repeat
        end

        coj = []; ndj = [];
        if wantmodes
          % ...
          
        end
        
      else, error(fprintf('unknown method %s', meth));
      end
    
      p.kj = kj; p.ej = ej; p.coj = coj;  % save internally (waste of memory?)
    end
    
    function A = fillfredholmop(p, k) % .... helper for solvespectrum
    % FILLFREDHOLMOP - fill l^2 normed (sqrtwei) layer-potential matrix at k
      p.setoverallwavenumber(k);
      A = p.fillbcmatrix;   % fill matrix mapping density values to field values
      N = size(A,1); if size(A,2)~=N, error('A must be square!'); end
      % rescale so l^2 norms represent L^2 norms, as Jim Bremer also does...
      A = A .* repmat(1./p.sqrtwei, [N 1]); % note sqrtwei already left-applied 
      A = -2*A;     % for Dirichlet, -2 turns (D-1/2) into (1-2D) for det
    end
      
    function [co nd e] = solvemodecoeffs(p, k, meth, o) % ........ single efunc
      if nargin<3, meth='fd'; end   % default method
      if nargin<4, o = []; end


    end
  end
  % --------------------------------------------------------------------
  methods(Static)    % these don't need EVP obj to exist to call them...
    %weylcount    - put in basic weyl routine
    %end
    smoothweylcount 
  end
end
