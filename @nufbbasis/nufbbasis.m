% NUFBBASIS - Create a basis of irregular Fourier-Bessel functions.
% This objects creates a basis of irregular Fourier-Bessel functions of
% degree nu*(1:N) (for Fourier-Bessel sine functions) or degree nu*(0:N)
% (for Fourier-Bessel cosine functions). A nufbbasis object is created by
%
% nufb=nufbbasis(origin,nu,offset,branch,N,k,opts)
% origin: The origin of the Fourier-Bessel functions
% nu    : The fractional order of the basis
% offset: The direction that corresponds to the angular variable theta=0
% branch: Direction of the branch cut. It must not point into the domain
% N     : Degree of the basis fct.
% opts  : Optional arguments. Currently supported is
%         - opts.type    : 's'  Create basis of Fourier-Bessel sine fct.
%                          'c'  Create basis of Fourier-Bessel cosine fct.
%                          'cs' Create basis of Fourier-Bessel cosine and
%                               sine functions 

classdef nufbbasis < handle & basis

    % Basis of irregular Fourier-Bessel functions.

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
        rescale_arg % argument (kR) to rescale the J bessels
    end

    methods
        function nufb = nufbbasis(origin,nu,offset,branch,N,k,opts)
            % 
            if nargin<7, opts=[]; end;
            if ~isfield(opts,'type'), opts.type='cs'; end;
            if nargin<6, k=NaN; end;
            if nargin<5, N=20; end;

            nufb.k=k;
            nufb.N=N;
            nufb.branch=branch;
            nufb.offset=offset;
            nufb.nu=nu;
            nufb.origin=origin;
            if opts.type=='s', nufb.type='s'; nufb.Nf=N; end;
            if opts.type=='c', nufb.type='c'; nufb.Nf=N+1; end;
            if strcmp(opts.type,'cs'), nufb.type='cs'; nufb.Nf=2*N+1; end;
            if ~isfield(opts,'rescale_rad'), opts.rescale_rad=0; end
            nufb.rescale_arg=k*opts.rescale_rad;

        end
        function [A, A1, A2] = eval(nufb,pts,opts)

            resc=(nufb.rescale_arg>0);
            % Evaluates the basis at a given set of points
            N=nufb.N; k=nufb.k; nu=nufb.nu; origin=nufb.origin;
            np=length(pts.x); % Number of points
            R=abs(pts.x-nufb.origin);
            ang=angle(-(pts.x-nufb.origin)./nufb.branch);
            offang=angle(-nufb.offset./nufb.branch);
            bes=besselj((0:N)*nu,k*R);
            if resc
               scfac=1./nufb.Jrescalefactors(0:N);
               bes=bes.*repmat(scfac,[numel(ang) 1]);
            end            
            if nufb.type=='s',
                %bes=besselj((1:N)*nu,k*R);
                s=sin(nu*(ang-offang)*(1:N));
                A=bes(:,2:end).*s;
            elseif nufb.type=='c',
                %bes=besselj((0:N)*nu,k*R);
                c=cos(nu*(ang-offang)*(0:N));
                A=bes.*c;
            elseif strcmp(nufb.type,'cs'),
                %bes=besselj((0:N)*nu,k*R);
                s=sin(nu*(ang-offang)*(1:N));
                c=cos(nu*(ang-offang)*(0:N));
                A=[bes.*c,bes(:,2:end).*s];
            end
            if nargout>1,
                if numel(find(R==0))>0,
                    warning('Computing x/y or normal derivatives of irregular Bessel functions at origin not implemented');
                end
                besr=k/2*(besselj(nu*(0:N)-1,k*R)-besselj(nu*(0:N)+1,k*R));
                if resc
                   scfac=1./nufb.Jrescalefactors(0:N);
                   besr=besr.*repmat(scfac,[numel(ang) 1]);
                end            
                if nufb.type=='s',
                    Ar=besr(:,2:end).*s;
                    At=nu*repmat(1:N,np,1).*bes(:,2:end).*cos(nu*(ang-offang)*(1:N));
                elseif nufb.type=='c',
                    Ar=besr.*c;
                    At=-nu*repmat(0:N,np,1).*bes.*sin(nu*(ang-offang)*(0:N));
                elseif strcmp(nufb.type,'cs'),
                    Ar=[besr.*c,besr(:,2:end).*s];
                    At=[-nu*repmat(0:N,np,1).*bes.*sin(nu*(ang-offang)*(0:N)),...
                        nu*repmat(1:N,np,1).*bes(:,2:end).*cos(nu*(ang-offang)*(1:N))];
                end
                ang0=angle(pts.x-origin); % Angle with respect to the original
                % coordinate system shifted by
                % origin
                cc=repmat(cos(ang0),1,nufb.Nf); ss=repmat(sin(ang0),1,nufb.Nf);
                RR=repmat(R,1,nufb.Nf);
                if nargout==2,
                    nx=repmat(real(pts.nx),1,nufb.Nf); ny=repmat(imag(pts.nx),1,nufb.Nf);
                    A1=Ar.*(nx.*cc+ny.*ss)+At.*(ny.*cc-nx.*ss)./RR;
                end
                if nargout==3,
                    A1=cc.*Ar-ss.*At./RR; A2=ss.*Ar+cc.*At./RR;
                end
            end
        end  % func


        function showgeom(b, opts) % .............. plotting of nufb basis geom
          if nargin<2, opts = []; end
          plot(real(nufb.origin), imag(nufb.origin), 'ro');
          str = sprintf('(nu=%.2g) ', b.nu);
          if isfield(opts, 'label'), str = [str opts.label]; end
          text(real(nufb.origin), imag(nufb.origin), str);
        end % func
        
        function sc = Jrescalefactors(nufb, n) % ... used to effect rescale_rad
        % JRESCALEFACTORS - given list of orders, return FB J-rescaling factors
          sc = abs(besselj(nufb.nu*n, min(nufb.nu*n, nufb.rescale_arg))); % note n can be list
          % the min here stops the J from getting to osc region (it stays 
          % at turning point)
        end


    end  % methods
end

            
            
            
            
            
            
