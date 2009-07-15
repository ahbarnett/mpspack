#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_errno.h>
#include <mex.h>
#include <matrix.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    int i;
    int m;
    int n;
    int nmin;
    int nmax;
    double* x;
    mxArray* besMatrix;
    double* besData;
    double* err;

    /* Check input and output data */


    if (nrhs!=3)
	mexErrMsgIdAndTxt("Utils:gslbesselj:nrhs",
			  "Wrong number of input arguments");

    if (!mxIsDouble(prhs[0]) ||
	mxGetNumberOfElements(prhs[1])!=1)
	mexErrMsgIdAndTxt("Utils:gslbesselj:notscalar",
			  "First input must be a scalar");

    if (!mxIsDouble(prhs[1]) ||
	mxGetNumberOfElements(prhs[1])!=1)
	mexErrMsgIdAndTxt("Utils:gslbesselj:notscalar",
			  "Second input must be a scalar");

    if (mxGetN(prhs[2])!=1)
	mexErrMsgIdAndTxt("Utils:gslbesselj:columnvector",
			  "Expected a column vector as input");

    /* Main routine */


    nmin=(int)mxGetScalar(prhs[0]);
    nmax=(int)mxGetScalar(prhs[1]);
    
    if (nmin<0)
	mexErrMsgIdAndTxt("Utils:gslbesselj:nmin",
			  "nmin must be a nonnegative integer");

    if (nmax<0)
	mexErrMsgIdAndTxt("Utils:gslbesselj:nmax",
			  "nmax must be a nonegative integer");

    if (nmax<nmin)
	mexErrMsgIdAndTxt("Utils:gslbesselj:nmax",
			  "nmax>=nmin required");


    x=mxGetPr(prhs[2]);
    m=mxGetM(prhs[2]);
    
    n=nmax-nmin+1;

    besMatrix=mxCreateDoubleMatrix(n,m,mxREAL);
    besData=mxGetPr(besMatrix);

    plhs[1]=mxCreateDoubleMatrix(m,1,mxREAL);
    err=mxGetPr(plhs[1]);

    gsl_set_error_handler_off();

    for (i=0; i<m; i++)
	err[i]=gsl_sf_bessel_Jn_array(nmin,nmax, x[i], &besData[i*n]);

    mexCallMATLAB(1,&plhs[0],1,&besMatrix,"transpose");

    return;

}
