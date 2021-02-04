#define _USE_MATH_DEFINES
#include <math.h>
#include "mex.h"

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *obs;
    double *p_start;
    double *p_trans;
    double *mu_emit;
    double *sigma_emit;
            
    double *outMatrix;
    
    mwSize dim;
    mwSize num_obs;
    mwSize num_states;

    /* sanity check inputs and extract the parameters to pass to the actual algorithm */
    
    /* check for proper number of arguments */
    
    if(nrhs!=5) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Inputs: obs [dim x #steps] start [1 x K] trans [K x K] mu [dim x K] sigma [1 x K]");
    }

    obs = mxGetPr(prhs[0]);
    p_start = mxGetPr(prhs[1]);
    p_trans = mxGetPr(prhs[2]);
    mu_emit = mxGetPr(prhs[3]);
    sigma_emit = mxGetPr(prhs[4]);
    
    dim = mxGetM(prhs[0]);
    num_obs = mxGetN(prhs[0]);
    num_states = mxGetN(prhs[1]);
    
    /* check input data type */
    
    if( !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3]) || !mxIsDouble(prhs[4])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","All inputs should be of type double.");
    }
    
    /* check that the dimensions of all matrices are consistent */
    
    if (dim != mxGetM(prhs[3])) {
        mexErrMsgIdAndTxt("MyToolbox:hmm:dim","obs and mu dimensions do not match.");
    }

    if (mxGetM(prhs[1]) != 1 || mxGetM(prhs[4]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:hmm:dim","p_start and sigma should be row vectors.");
    }

    if (mxGetN(prhs[1]) != num_states || mxGetM(prhs[2]) != num_states || 
            mxGetN(prhs[2]) != num_states || mxGetN(prhs[3]) != num_states || mxGetN(prhs[4]) != num_states) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","K state dimensions of p_start, p_trans, mu, and sigma are not consistent");
    }

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);

    /* get a pointer to the real data in the output matrix */
    outMatrix = mxGetPr(plhs[0]);

    /* run the routine */
    hmm_forward(obs, p_start, p_trans, mu_emit, sigma_emit, 
        dim, num_obs, num_states, outMatrix);
}
