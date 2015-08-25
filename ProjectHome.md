MPSpack is a user-friendly and fully object-oriented MATLAB toolbox that implements the method of particular solutions, nonpolynomial FEM, the method of fundamental solutions, and integral equation methods, for the efficient and highly accurate solution of Laplace eigenvalue problems, interior/exterior Helmholtz boundary-value problems (e.g. wave scattering), periodic diffraction problems, and related PDE problems, on piecewise-homogeneous 2D domains.

We released version 1.0 in September 2009. As of January 2014 the current version 1.33 includes Kress corner parametrization for layer potentials, and finally some documentation and worked examples for eigenvalue problems, including the Neumann case.

Please see the **Downloads** page for archives of the package, the Manual
which has installation instructions, and the all-important Tutorial.
See the **Source** page for how to download via `svn` (subversion); this gives access to the current development version.

This material is based upon work supported by the National Science Foundation under grants
DMS-0811005 and DMS-1216656, and Engineering and Physical Sciences Research Council Grant EP/H00409/1, among others. Since 2009 it has been supported by Alex Barnett alone. It includes codes by V. Rokhlin, A. Pataki, Z. Gimbutas, S. Hawkins, and several others.

_First is an example showing scattering from a square, accurate to 10 digits,
computed in a few seconds on a laptop.
Spectral convergence is achieved using the following ingredients: decomposition into subdomains (nonpolynomial FEM), fractional-order Fourier-Bessel expansions around corner singularities, and an exterior fundamental solutions representation. With MPSpack this
needs no more than 20 lines of code._

![http://math.dartmouth.edu/~ahb/images/sqscatt2_cut_sm.png](http://math.dartmouth.edu/~ahb/images/sqscatt2_cut_sm.png)

_Next is an example showing the first 45 Dirichlet eigenmodes of a smooth planar domain, computed to 12 digit accuracy and evaluated on a grid of 3600 points, in around 1 second per mode. Convergence is again spectral, using layer potential representations, and analytic root-finding on a Fredholm determinant. With MPSpack the code is 9 lines long._

![http://math.dartmouth.edu/~ahb/images/rf_45modes.png](http://math.dartmouth.edu/~ahb/images/rf_45modes.png)

_Finally we show an entertaining example of acoustic (sound-hard) scattering from some digitized letter shapes, computed recently to 10 digit accuracy by Perrin Meyer.
The wave is incident from about 4 o'clock:_
![http://math.dartmouth.edu/~ahb/images/hny2014_perrin.png](http://math.dartmouth.edu/~ahb/images/hny2014_perrin.png)