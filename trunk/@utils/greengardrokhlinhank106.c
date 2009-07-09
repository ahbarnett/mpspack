/* MEX interface to Greengard-Rokhlin precomputed Chebychev lookup for Hankels.

   Based on greengardrokhlinhank103.c
   Must be linked to hank103.o and hank106.o; see Makefile.

   barnett 9/5/08
 */

#include <mex.h>
#include <matrix.h>

extern void hank103_(double *a, double *b, double *c, int *d);
extern void hank106_(double *a, double *b, double *c, int *d);
extern void hank106datagen_(double *rk, int *ier);

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])
{
  int n, m, j;
  double *zr, *zi, *H0r, *H0i, *H1r, *H1i;  /* ptrs to real/imag array parts */
  double z[2], H0[2], H1[2];                /* non-array I/O registers */
  double rk[2] = {1.0,0.0};                 /* fake wavenumber for datagen */
  int ier, ifexpon = 1; /* if 1 unscaled; 0 hank103 scales by  e^{-i \cdot z} */

  if (nrhs != 1)
    mexErrMsgTxt("greengardrokhlinhank106: must have one input argument");
  if (nlhs > 2)
    mexErrMsgTxt("greengardrokhlinhank106: must have no more than two output arguments");
  if (mxIsClass(prhs[0],"sparse") || mxIsChar(prhs[0]))
    mexErrMsgTxt("greengardrokhlinhank106: input must be full and nonstring");

  zr = mxGetPr(prhs[0]); zi = mxGetPi(prhs[0]);     /* pointer to input (RHS) */

  /* allocate output (LHS) arrays and H?? as pointers to them... */
  m = mxGetM(prhs[0]); n = mxGetN(prhs[0]);
  plhs[0] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
  plhs[1] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
  H0r = mxGetPr(plhs[0]); H0i = mxGetPi(plhs[0]);
  H1r = mxGetPr(plhs[1]); H1i = mxGetPi(plhs[1]);

  hank106datagen_(rk, &ier);    /* generate lookup table (in static arrays) */

  if (!mxIsComplex(prhs[0]))            /* input is real, has no imag part */
    for (j=0;j<n*m;++j) {               /* note looping over array in one go */
      z[0] = zr[j]; z[1] = 0.0;         /* copy real part to register, imag=0 */
      hank106_(z, H0, H1, &ifexpon);    /* pass by reference to fortran */
      H0r[j] = H0[0]; H0i[j] = H0[1]; H1r[j] = H1[0]; H1i[j] = H1[1];
    }
  else                                  /* input is complex */
    for (j=0;j<n*m;++j) {
      z[0] = zr[j]; z[1] = zi[j];       /* copy to register so contiguous */
      hank106_(z, H0, H1, &ifexpon);
      H0r[j] = H0[0]; H0i[j] = H0[1]; H1r[j] = H1[0]; H1i[j] = H1[1];
    }

  return;
}
