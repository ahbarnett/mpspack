% MPSpack: Method of Particular Solutions Toolbox.
%          A MATLAB toolbox to solve Helmholtz problems with particular and
%          fundamental solution methods
%
%          authors: Alex Barnett and Timo Betcke
%
%   This is a summary of commands; see help command for more detail on each
%
% Pointsets.
%   pointset    - create a pointset
%   plot        - draw a pointset
%
% Segments.
%   segment     - create segment object with quadrature pts, weights, normals
%   polyseglist - create closed segment object list from CCW polygon vertices
%   plot        - draw a segment
%   showsegments - draw connected segment list with signs
%   testsegment - test program for the segment class
%
% Domains.
%   domain      - create interior/exterior domain from segment and sign lists
%   plot        - draw a domain, with options for various geometric aspects
%   testdomain  - test program for the domain class
%
% Basis sets.
%   regfbbasis  - create a regular Fourier-Bessel basis set object
%   eval        - evaluate any basis set object on a point set (eg segment)
