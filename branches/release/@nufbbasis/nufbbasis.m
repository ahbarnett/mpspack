classdef nufbbasis < handle & basis

% NUFBBASIS - create a fractional-order Fourier-Bessel basis set (corner exp)
%
% b = NUFBBASIS(origin, nu, offset, branch, N) creates a basis of fractional-
%   order Fourier-Bessel functions appropriate for expansion of the Helmholtz
%   equation in a wedge of angle pi/nu > 0. The orders are nu*(1:N) (for sine
%   angular functions) or nu*(0:N) (for cosine angular functions). The other
%   arguments are, 
%     origin: The origin of the Fourier-Bessel functions (as complex number)
%     offset: The direction that corresponds to the angular variable theta=0
%     branch: Direction of the branch cut (point on the unit circle);
%             it must not point into the affected domain.
%     N     : degree of the basis set
%
%   The wedge lies counterclockwise from offset, spanning pi/nu in angle.
%   As with all basis types, the wavenumber k is determined by that of the
%   affected domain. k=0 gives fractional-power Laplace equation solutions.
%
% b = NUFBBASIS(origin, nu, offset, branch, N, opts) as above but permits
%   optional arguments to be selected. Currently supported is
%         - opts.type    : 's'  Create basis of Fourier-Bessel sine fct.
%                          'c'  Create basis of Fourier-Bessel cosine fct.
%                          'cs' Create basis of Fourier-Bessel cosine and
%                               sine functions (default)
%   The reason for having 'cs' rather than using separate 's' and 'c' objects
%   is a factor of 2 in speed: the set of Bessel evaluations is reused.
%
% Issues/notes:
%  * need mixed-type D-N wedge angular functions too.
%  * add switch to GSL bessel_nu
%  * for slit domains nu=1/2, user needs to jiggle the branch cut either side?
%
% Also see: DOMAIN.ADDNUFBBASIS

% Copyright (C) 2008, 2009, Timo Betcke, Alex Barnett

    properties
        origin   % Origin of the Fourier-Bessel fct.
        nu       % Real parameter that determines the fractional order of
                 % the Bessel fct.
        branch   % Complex angle that determines the branch cut
        offset   % Complex angle that determines the zero line of the Bessel
                 % fct.
        type     % 'c': Fourier-Bessel cosine basis
                 % 's': Fourier-Bessel sine basis
                 % 'cs': Fourier-Bessel cosine/sine basis (default)
        rescale_rad % radius used to rescale the J bessels or power laws
    end

    methods
        function nufb = nufbbasis(origin,nu,offset,branch,N,opts)
            if nargin<5, N=20; end;
            if nargin<6, opts=[]; end;
            if ~isfield(opts,'type'), opts.type='cs'; end;
            if ~isfield(opts,'rescale_rad'), opts.rescale_rad=0; end
            if ~isfield(opts,'nmultiplier'), opts.nmultiplier=1; end

            nufb.nmultiplier=opts.nmultiplier; nufb.updateN(N);
            nufb.branch=branch; nufb.offset=offset;
            nufb.nu=nu; nufb.origin=origin;
            nufb.type=opts.type; nufb.rescale_rad=opts.rescale_rad;
            nufb.Nf;                 % raises an error if unknown type
        end
        
        function Nf = Nf(b) % ......................... # degrees of freedom
          switch b.type
           case 's'
            Nf = b.N;
           case 'c'
            Nf = b.N + 1;
           case 'cs'
            Nf = 2*b.N + 1;
           otherwise
            error('unknown cos/sine angular type in nufbbasis!');
          end
        end

        function [A, A1, A2] = eval(nufb,pts,opts)
        % EVAL - evaluates fractional Fourier-Bessel basis at a set of points
        %
        % A = EVAL(p) where p is a pointset object containing M points, returns
        %   a M-by-Nf matrix whose ij'th entry is Phi_j(z_i), where Phi_j is
        %   the jth basis function, and z_i the ith point. Nf is the number
        %   of degrees of freedom in the basis set object.
        %        
        % [A An] = EVAL(p) also returns matrix An whose ij'th entry is
        %   d/dn_i Phi_j(z_i), the derivative in the ith normal direction in
        %   the pointset.
        %
        % [A Ax Ay] = EVAL(p) returns A as above, and matrices of basis function
        %    derivatives in the x- and y-directions. That is, Ax has ij'th
        %    entry d/dx Phi_j(z_i) while Ay has ij'th entry d/dy Phi_j(z_i)
        %
        % Also see: POINTSET, NUFBBASIS

            resc=(nufb.rescale_rad>0); % whether to rescale or not
            N = nufb.N;
            k = nufb.k;                 % NB this is now a method not property
            if isnan(k), error('nufbbasis: no domain wavenumber defined!'); end
            nu=nufb.nu; origin=nufb.origin;
            np=length(pts.x);                        % Number of points
            R=abs(pts.x-origin);
            ang=angle(-(pts.x-origin)./nufb.branch); % col vec of angles
            offang=angle(-nufb.offset./nufb.branch);
            % do the expensive evaluation...
            if k==0, bes = repmat(R, [1 N+1]).^repmat((0:N)*nu, [np 1]);
            else, bes=besselj((0:N)*nu,k*R);
            end
            if resc                                    % rescale by orders
              scfac = 1./nufb.Jrescalefactors(0:N);    % uses rescale_arg
              bes=bes.*repmat(scfac,[np 1]);
            end
            if nufb.type=='s',
                s=sin(nu*(ang-offang)*(1:N));  % note * is outer product
                A=bes(:,2:end).*s;
            elseif nufb.type=='c',
                c=cos(nu*(ang-offang)*(0:N));
                A=bes.*c;
            elseif strcmp(nufb.type,'cs'),
                s=sin(nu*(ang-offang)*(1:N));
                c=cos(nu*(ang-offang)*(0:N));
                A=[bes.*c,bes(:,2:end).*s];
            end
            
            if nargout>1       % --------------- derivs wanted ----------
              i = find(R==0); R(i) = 1;     % zero-R indices, use dummy R vals
              if k==0
                besr = repmat((0:N)*nu, [np 1]).* ...
                       repmat(R, [1 N+1]).^repmat((0:N)*nu-1, [np 1]);;
              else  % fix this, which involves 2 new besselnu calls...
                besr=k/2*(besselj(nu*(0:N)-1,k*R)-besselj(nu*(0:N)+1,k*R));
              end
              if resc, besr=besr.*repmat(scfac,[np 1]); end  % rescale
              if nufb.type=='s'
                nus = (1:N)*nu;        % list of orders
                Ar=besr(:,2:end).*s;
                At=repmat(nus,np,1).*bes(:,2:end).*cos(nu*(ang-offang)*(1:N));
              elseif nufb.type=='c'
                nus = (0:N)*nu;
                Ar=besr.*c;
                At=-repmat(nus,np,1).*bes.*sin(nu*(ang-offang)*(0:N));
              elseif strcmp(nufb.type,'cs'),
                nus = [(0:N)*nu (1:N)*nu];
                Ar=[besr.*c,besr(:,2:end).*s];
                % edited to reuse c and s arrays from above...
                At=[-repmat((0:N)*nu,np,1).*bes.*[0*R s],...
                    repmat((1:N)*nu,np,1).*bes(:,2:end).*c(:,2:end)];
              end
              ang0=angle(pts.x-origin); % Angle with respect to the original
                                        % coordinate system shifted by
                                        % origin (i.e. without branch rotation)
              cc=repmat(cos(ang0),1,nufb.Nf); ss=repmat(sin(ang0),1,nufb.Nf);
              invRR=repmat(1./R,1,nufb.Nf); % since dividing is slow
              nn = find(nus==1);        % indices of orders that are 1

              if nargout==2      % normal directional deriv ------
                nx=repmat(real(pts.nx),1,nufb.Nf);
                ny=repmat(imag(pts.nx),1,nufb.Nf);
                A1=Ar.*(nx.*cc+ny.*ss)+At.*(ny.*cc-nx.*ss).*invRR;
                if ~isempty(i)  % overwrite the zero-R derivs... (pain for nu=1)
                  nnor = repmat(pts.nx./nufb.offset, size(nn)); % nang vs off
                  A1(i,find(nus>1 | nus==0)) = 0;    % order 0 is special
                  A1(i,find(nus<1 & nus~=0)) = Inf;  % divergence of all derivs
                  if nufb.type=='s'
                    A1(i,nn) = imag(nnor(i,:));            % nu=1 sin term
                  elseif nufb.type=='c'
                    A1(i,nn) = real(nnor(i,:));            % nu=1 cos term
                  elseif strcmp(nufb.type,'cs')
                    cnn = find(nus==1 & (1:nufb.Nf)<=N+1);
                    snn = find(nus==1 & (1:nufb.Nf)>N+1);
                    A1(i,cnn) = real(repmat(pts.nx(i)./nufb.offset, size(cnn)));
                    A1(i,snn) = imag(repmat(pts.nx(i)./nufb.offset, size(cnn)));
                  end
                  if k>0, A1(i,nn) = k/2 * A1(i,nn); end
                  if ~isempty(nn) & resc
                    A1(i,nn) = A1(i,nn) / nufb.Jrescalefactors(1); end
                end
                
              elseif nargout==3     % x- and y-derivs --------
                A1=cc.*Ar-ss.*At.*invRR; A2=ss.*Ar+cc.*At.*invRR;
                if ~isempty(i)  % overwrite the zero-R derivs... (pain for nu=1)
                  A1(i,find(nus>1 | nus==0)) = 0;
                  A2(i,find(nus>1 | nus==0)) = 0;
                  A1(i,find(nus<1 & nus~=0)) = Inf;
                  A2(i,find(nus<1 & nus~=0)) = Inf;
                  if nufb.type=='s'
                    A1(i,nn) = -imag(nufb.offset);            % nu=1 sin term
                    A2(i,nn) = real(nufb.offset);
                  elseif nufb.type=='c'
                    A1(i,nn) = real(nufb.offset);            % nu=1 cos term
                    A2(i,nn) = imag(nufb.offset);
                  elseif strcmp(nufb.type,'cs')
                    cnn = find(nus==1 & (1:nufb.Nf)<=N+1);
                    snn = find(nus==1 & (1:nufb.Nf)>N+1);
                    A1(i,cnn) = real(nufb.offset);
                    A2(i,cnn) = imag(nufb.offset);
                    A1(i,snn) = -imag(nufb.offset);
                    A2(i,snn) = real(nufb.offset);
                  end
                  if k>0, A1(i,nn) = k/2 * A1(i,nn); A2(i,nn) = k/2 * A2(i,nn);
                  end
                  if ~isempty(nn) & resc, scfaco = 1./nufb.Jrescalefactors(1);
                    A1(i,nn) = A1(i,nn) * scfaco; A2(i,nn) = A2(i,nn) * scfaco;
                  end                 
                end
              end
            end
        end  % func

        function h = showgeom(nufb, opts) % ....... plotting of nufb basis geom
        % SHOWGEOM - show corner wedge location, geometry and branch cut
          if nargin<2, opts = []; end
          x = real(nufb.origin); y = imag(nufb.origin); h = plot(x,y,'ro');
          str = sprintf('(nu=%.2g) ', nufb.nu);
          if isfield(opts, 'label'), str = [str opts.label]; end
          h = [h; text(x,y, str)]; hold on;
          angfrac = 0:.05:1;            % see code in domain.plot corners:
          t = nufb.offset * exp(1i*pi/nufb.nu.*angfrac); % angles on unit circle
          l = 0.3; % filled polygon patch...
          h = [h; patch([x+l*real(t) x], [y+l*imag(t) y], 'r')];
          h = [h; plot([x+l*real(nufb.branch) x], [y+l*imag(nufb.branch) y],...
                       'r--', 'linewidth', 5)];
        end % func
        
        function sc = Jrescalefactors(nufb, n) % ... used to effect rescale_rad
        % JRESCALEFACTORS - given list of orders, return FB J-rescaling factors
        %
        %  Also now works for k=0, ie harmonic polynomials.
          k = nufb.k;               % method gets wavenumber
          if k==0                   % polynomial case
            sc = nufb.rescale_rad.^n;
          else                      % Bessel case
            sc = abs(besselj(nufb.nu*n, min(nufb.nu*n, k*nufb.rescale_rad)));
            % Note n can be list.
            % The min here stops the J from getting to osc region (it stays 
            % at turning point)
          end
        end


    end  % methods
end

            
            
            
            
            
            
