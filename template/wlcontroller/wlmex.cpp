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
	Eigen::Vector4d uoute = Eigen::Map<Eigen::Vector4d>(uout);
	// init
	wlc.setLimits(u_t(90., -0.5, -0.2, -0.1), u_t(160., 0.5, 0.2, 0.1), u_t(5., 0.01, 0.01, 0.01));
	wlc.setWeight(MC6(Qdiag).cast<float>());
	uoute = wlc.update(MC4(u0).cast<float>(), MC6(p0).cast<float>(), MC6(h0).cast<float>(), MC6(pdes).cast<float>(), MC6(kpmom).cast<float>()).cast<double>();
}

/* The gateway function */
extern "C"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double multiplier;              /* input scalar */
	double *inMatrix;               /* 1xN input matrix */
	size_t ncols;                   /* size of matrix */
	double *outMatrix;              /* output matrix */

	/* check for proper number of arguments */
	if(nrhs!=2) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
	}
	if(nlhs!=1) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
	}
	/* make sure the first input argument is scalar */
	if( !mxIsDouble(prhs[0]) || 
			mxIsComplex(prhs[0]) ||
			mxGetNumberOfElements(prhs[0])!=1 ) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Input multiplier must be a scalar.");
	}

	/* make sure the second input argument is type double */
	if( !mxIsDouble(prhs[1]) || 
			mxIsComplex(prhs[1])) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
	}

	/* check that number of rows in second input argument is 1 */
	if(mxGetM(prhs[1])!=1) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
	}

	/* get the value of the scalar input  */
	multiplier = mxGetScalar(prhs[0]);

	/* create a pointer to the real data in the input matrix  */
	inMatrix = mxGetPr(prhs[1]);

	/* get dimensions of the input matrix */
	ncols = mxGetN(prhs[1]);

	/* create the output matrix */
	plhs[0] = mxCreateDoubleMatrix(1,(mwSize)ncols,mxREAL);

	/* get a pointer to the real data in the output matrix */
	outMatrix = mxGetPr(plhs[0]);

	/* call the computational routine */
	// arrayProduct(multiplier,inMatrix,outMatrix,(mwSize)ncols);
}
