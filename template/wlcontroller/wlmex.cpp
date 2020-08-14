/**
 * @file wlmex.cpp
 * @author Avik De (avikde@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-08-14
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include <mex.h>
#include "wlcontroller.hpp"

// Global objects are OK https://stackoverflow.com/questions/39756108/matlab-mex-files-with-global-data-in-c
WLController wlc;
typedef Eigen::Map<const Eigen::Matrix<double, 6, 1> > MC6;
typedef Eigen::Map<const Eigen::Matrix<double, 4, 1> > MC4;

static void wlcWrapper(double uout[/* 4 */], const double u0[/* 4 */], const double p0[/* 6 */], const double h0[/* 6 */], const double pdes[/* 6 */], const double kpmom[/* 6 */], const double Qdiag[/* 6 */]) {
	// init
	wlc.setLimits(u_t(90., -0.5, -0.2, -0.1), u_t(160., 0.5, 0.2, 0.1), u_t(5., 0.01, 0.01, 0.01));
	wlc.setWeight(MC6(Qdiag).cast<float>());

	// Actual solve
	Eigen::Vector4d uoute = wlc.update(MC4(u0).cast<float>(), MC6(p0).cast<float>(), MC6(h0).cast<float>(), MC6(pdes).cast<float>(), MC6(kpmom).cast<float>()).cast<double>();

	// For some reason Eigen::Map to assign didn't work here
	memcpy(uout, uoute.data(), 4 * sizeof(double));
}

// Entry point
extern "C"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// /* check for proper number of arguments */
	if(nrhs < 6) {
		mexErrMsgIdAndTxt("wlmex:nrhs","Six inputs required.");
	}
	if(nlhs != 1) {
		mexErrMsgIdAndTxt("wlmex:nlhs","One output required.");
	}
	// /* make sure the first input argument is scalar */
	if (!mxIsDouble(prhs[0]) || mxGetNumberOfElements(prhs[0]) != 4) {
		mexErrMsgIdAndTxt("wlmex", "u0 must have 4 doubles");
	}
	if (!mxIsDouble(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 6) {
		mexErrMsgIdAndTxt("wlmex", "p0 must have 6 doubles");
	}
	if (!mxIsDouble(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 6) {
		mexErrMsgIdAndTxt("wlmex", "h0 must have 6 doubles");
	}
	if (!mxIsDouble(prhs[3]) || mxGetNumberOfElements(prhs[3]) != 6) {
		mexErrMsgIdAndTxt("wlmex", "pdes must have 6 doubles");
	}
	if (!mxIsDouble(prhs[4]) || mxGetNumberOfElements(prhs[4]) != 6) {
		mexErrMsgIdAndTxt("wlmex", "kpmom must have 6 doubles");
	}
	if (!mxIsDouble(prhs[5]) || mxGetNumberOfElements(prhs[5]) != 6) {
		mexErrMsgIdAndTxt("wlmex", "Qdiag must have 6 doubles");
	}
	
	// Create room for output
	plhs[0] = mxCreateDoubleMatrix(1, 4, mxREAL);

	wlcWrapper(mxGetPr(plhs[0]), mxGetPr(prhs[0]), mxGetPr(prhs[1]), mxGetPr(prhs[2]), mxGetPr(prhs[3]), mxGetPr(prhs[4]), mxGetPr(prhs[5]));
}
