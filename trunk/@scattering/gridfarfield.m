% GRIDFARFIELD - setup a grid and compute the far field on the grid.
%
% u = p.gridfarfield(opts) computes the far field u for the scattering problem
% p. The grid spacing is specified by opts.dx.
%
% [u, theta] = p.gridfarfield(opts) also returns the grid points.
%
% Copyright (C) 2014 Stuart C. Hawkins

function [u theta] = gridfarfield(self,opts)

%-----------------------------------
% setup the grid... this is a set
% of points on the unit circle parametrised
% by theta
%-----------------------------------

% set theta points
theta = linspace(0,2*pi,ceil(2*pi/opts.dx));

% get complex values points on the unit circle
z = exp(1i*theta);

%-----------------------------------
% compute the far field
%-----------------------------------

u = self.pointfarfield(z);