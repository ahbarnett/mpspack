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

    nmin=(int)mxGetScalar(prhs[0]);
    nmax=(int)mxGetScalar(prhs[1]);
    
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
