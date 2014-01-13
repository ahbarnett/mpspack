% SHOWFARFIELD - plot the farfield in a figure.
%
% showfarfield(pr) plots the absolute value of the farfield of the solution
% of the scattering problem pr.
%
% showfarfield(pr,opts) specifies an options structure with optional
% fields:
%   opts.real = true, plots the real part of the farfield
%   opts.imag = true, plots the imaginary part of the farfield
%   opts.abs = true, plots the absolute value of the farfield
%   opts.dB = true, plots the intensity of the farfield in decibels (engineers)
%
% [u,theta] = showfarfield(pr,opts) also returns the farfield values u and 
% the corresponding azimuth angles theta.
%
% Copyright (C) 2014 Stuart C. Hawkins, tweaks by Alex Barnett

function [u theta] = showfarfield(self,opts)

%-----------------------------------
% extract parameters
%-----------------------------------

% set default for opts
if nargin < 2
    opts = [];
end

% default line style
if ~isfield(opts,'plotopts')
    opts.plotopts = 'b-';
end

% default: absolute value not computed
if ~isfield(opts,'abs')
    opts.abs = 0;
end

% default: imaginary part not computed
if ~isfield(opts,'imag')
    opts.imag = 0;
end

% default: real part not computed
if ~isfield(opts,'real')
    opts.real = 0;
end

% default: not dB
if ~isfield(opts,'dB')
    opts.dB = 0;
end

% now... set default behaviour for the plot
if ~opts.real & ~opts.imag
    opts.abs = 1;
end

%-----------------------------------
% get the far field
%-----------------------------------

[u,theta] = self.gridfarfield(opts);

%-----------------------------------
% plot the data
%-----------------------------------

% get the hold state
hold_state = ishold;

if opts.dB      % added Barnett
  plot(theta,10*log10(abs(u).^2), opts.plotopts);
  ylabel('dB'); hold on

elseif opts.abs
  plot(theta,abs(u),opts.plotopts)
  hold on
end

if opts.real & opts.imag
    plot(theta,real(u),opts.plotopts,theta,imag(u),opts.plotopts)
    hold on
end

if opts.real
    plot(theta,real(u),opts.plotopts)
    hold on
end

if opts.imag
    plot(theta,imag(u),opts.plotopts)
    hold on
end    

if ~hold_state
    hold off
end
