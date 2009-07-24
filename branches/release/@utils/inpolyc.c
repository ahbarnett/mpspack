/* MEX interface to inpolygon implementation code, MPSpack toolbox.

Copyright (C) 2008, 2009, Timo Betcke, Alex Barnett
*/

#include <mex.h>

/* Original C-code acknowledgments (pnpoly):

Copyright (c) 1970-2003, Wm. Randolph Franklin

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

   1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimers.
   2. Redistributions in binary form must reproduce the above
   copyright notice in the documentation and/or other materials
   provided with the distribution.  3. The name of W. Randolph
   Franklin may not be used to endorse or promote products derived
   from this Software without specific prior written permission.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

int pnpoly(int nvert, double *vertx, double *verty, double testx, double testy)
/* Argument	Meaning
   nvert 	Number of vertices in the polygon. Whether to repeat the
                first vertex at the end is discussed below.
   vertx, verty Arrays containing the x- and y-coordinates of the polygon's
                vertices.
   testx, testy	X- and y-coordinate of the test point.

   This is the core computational routine acting on a single test point.
*/
{
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *px, *py, *vx, *vy;
  mwSize N, nv;
  int j;
  int *i;   /* still not convinced this is correct, but ok */

  if (nrhs!=2) mexErrMsgTxt("inpolyc.c: Wrong number of input arguments");
  if (nlhs>=2) mexErrMsgTxt("inpolyc.c: Wrong number of output arguments");

  /* handle the empty p input case */
  if (mxGetN(prhs[0])==0) { 
    plhs[0] = mxCreateNumericMatrix(0, 0, mxINT32_CLASS, 0);
    return;
  }

  if (mxGetN(prhs[0])!=1)
    mexErrMsgTxt("inpolyc.c: p must be column vector");
  N = (int)mxGetM(prhs[0]);    /* number of test points */
  px = mxGetPr(prhs[0]);
  if (mxIsComplex(prhs[0]))
    py = mxGetPi(prhs[0]);
  else                     /* create zero imag part */
    py = (double *) mxCreateNumericMatrix(N, 1, mxDOUBLE_CLASS, 0);

  if (!mxIsEmpty(prhs[1]) && mxGetN(prhs[1])!=1)
    mexErrMsgTxt("inpolyc.c: v must be column vector (or empty)");
  nv = (int)mxGetM(prhs[1]);    /* number of vertices */
  vx = mxGetPr(prhs[1]);
  if (mxIsComplex(prhs[1]))
    vy = mxGetPi(prhs[1]);
  else                     /* create zero imag part */
    vy = (double *) mxCreateNumericMatrix(nv, 1, mxDOUBLE_CLASS, 0);

  /* allocate output... column vector of 4-byte integers */
  plhs[0] = mxCreateNumericMatrix(N, 1, mxINT32_CLASS, 0);
  i = (int *) mxGetPr(plhs[0]);
  
  for (j=0; j<N; ++j)        /* repeated calls to check each pts one by one */
    i[j] = pnpoly(nv, vx, vy, px[j], py[j]);

  /* printf("%d %d\n", N,nv); */
  return;
}

