/**
 * @file wlmexc.c
 * @author Avik De (avikde@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-08-14
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include <mex.h>
#include "wlcontroller.h"

// Entry point
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// /* check for proper number of arguments */
	if(nrhs < 3) {
		mexErrMsgIdAndTxt("wlControllerUpdate:nrhs","Six inputs required.");
	}
	if(nlhs != 1) {
		mexErrMsgIdAndTxt("wlControllerUpdate:nlhs","One output required.");
	}
	// /* make sure the first input argument is scalar */
	if (!mxIsSingle(prhs[0]) || mxGetNumberOfElements(prhs[0]) != 4) {
		mexErrMsgIdAndTxt("wlControllerUpdate:notSingle", "u0 must have 4 floats");
	}
	if (!mxIsSingle(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 6) {
		mexErrMsgIdAndTxt("wlControllerUpdate:notSingle", "h0 must have 6 floats");
	}
	if (!mxIsSingle(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 6) {
		mexErrMsgIdAndTxt("wlControllerUpdate:notSingle", "pddes must have 6 floats");
	}
	
	// Create room for output
	plhs[0] = mxCreateNumericMatrix(1, 4, mxSINGLE_CLASS, mxREAL);

	wlControllerUpdate(mxGetPr(plhs[0]), mxGetPr(prhs[0]), mxGetPr(prhs[1]), mxGetPr(prhs[2]));
}
