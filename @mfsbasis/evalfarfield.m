% EVALFARFIELD - evaluate the far field from point sources.
%
% A = b.evalfarfield(pts) evaluate the far field of the mfsbasis b at
% points specified by pts.x.
%
% Copyright (C) 2014 Stuart C. Hawkins

function A = evalfarfield(self,pts)

%-----------------------------------
% extract parameters
%-----------------------------------

% number of sources
N = self.N;

% number of points
M = numel(pts.x);

% wavenumber
k = self.k;

% coupling parameter
eta = self.eta;

%-----------------------------------
% setup
%-----------------------------------

% get normal at the y points
ny = repmat(self.ny, [M 1]);

% get x points
tx = pts.x(:);
x = repmat(tx, [1 N]);

% get xh points on unit circle from x points
xh=x./abs(x);

% get y points
ty = self.y(:).';
y = repmat(ty, [M 1]);

%-----------------------------------
% compute matrix A with A(i,j) being
% the far field at x(j) due to source
% at x(i)
%-----------------------------------

if isinf(eta)
    
    A=(1i/4)*(1/sqrt(pi*k))*(1-1i)*exp(-1i*k*real(xh.*conj(y)));


elseif eta==0

    A=k*(1i/4)*(1/sqrt(pi*k))*(-1-1i)*real(xh.*conj(ny)).*exp(-1i*k*real(xh.*conj(y)));
    
else
    
    A=-eta*(1/4)*(1/sqrt(pi*k))*(1-1i)*exp(-1i*k*real(xh.*conj(y))) + ...
        k*(1i/4)*(1/sqrt(pi*k))*(-1-1i)*real(xh.*conj(ny)).*exp(-1i*k*real(xh.*conj(y)));
    
end
