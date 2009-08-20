#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_errno.h>
#include <mex.h>
#include <matrix.h>
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    int i;
    int vsize;
    int xsize;
    
    double* v;
    double* x;
    double* err;
    double* besData;
    double* pr;

    mxArray* rhs[3];
    

    /* Check input and output data */


/*     if (nrhs!=3) */
/* 	mexErrMsgIdAndTxt("Utils:gslbesselj:nrhs", */
/* 			  "Wrong number of input arguments"); */

/*     if (!mxIsDouble(prhs[0]) || */
/* 	mxGetNumberOfElements(prhs[1])!=1) */
/* 	mexErrMsgIdAndTxt("Utils:gslbesselj:notscalar", */
/* 			  "First input must be a scalar"); */

/*     if (!mxIsDouble(prhs[1]) || */
/* 	mxGetNumberOfElements(prhs[1])!=1) */
/* 	mexErrMsgIdAndTxt("Utils:gslbesselj:notscalar", */
/* 			  "Second input must be a scalar"); */

/*     if (mxGetN(prhs[2])!=1) */
/* 	mexErrMsgIdAndTxt("Utils:gslbesselj:columnvector", */
/* 			  "Expected a column vector as input"); */

    /* Main routine */


    vsize=(int)mxGetN(prhs[0]);
    xsize=(int)mxGetM(prhs[1]);
    
/*     if (nmin<0) */
/* 	mexErrMsgIdAndTxt("Utils:gslbesselj:nmin", */
/* 			  "nmin must be a nonnegative integer"); */

/*     if (nmax<0) */
/* 	mexErrMsgIdAndTxt("Utils:gslbesselj:nmax", */
/* 			  "nmax must be a nonegative integer"); */

/*     if (nmax<nmin) */
/* 	mexErrMsgIdAndTxt("Utils:gslbesselj:nmax", */
/* 			  "nmax>=nmin required"); */


    v=mxGetPr(prhs[0]);
    x=mxGetPr(prhs[1]);

    /* Fill up besData with copies of x */

/*     for (i=0; i<vsize; i++) */
/* 	memcpy(&besData[i*xsize],x,xsize*sizeof(double)); */

    rhs[0]=prhs[1];
    rhs[1]=mxCreateDoubleMatrix(1,1,mxREAL);
    rhs[2]=mxCreateDoubleMatrix(1,1,mxREAL);

    pr=mxGetPr(rhs[1]);
    *pr=1;
    
    pr=mxGetPr(rhs[2]);
    *pr=vsize;

    mexCallMATLAB(1,&plhs[0],3,rhs,"repmat");

    besData=mxGetPr(plhs[0]);
    
    /* Create Matrix for error handler */

    plhs[1]=mxCreateDoubleMatrix(1,vsize,mxREAL);
    err=mxGetPr(plhs[1]);

    gsl_set_error_handler_off();

    for (i=0; i<vsize; i++)
	err[i]=gsl_sf_bessel_sequence_Jnu_e(v[i],GSL_PREC_DOUBLE,xsize,&besData[i*xsize]);

    return;

}
