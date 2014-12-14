#include "mex.h"
#include "rodrigues.cpp"

double *R, *w;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	if(nrhs!=1) mexErrMsgTxt("Needs 1 input argument");
// 	if(nlhs!=1) mexErrMsgTxt("Must have 1 output variable");
	
	plhs[0] = mxCreateDoubleMatrix(3,1,mxREAL);
	w = (double*)mxGetPr(plhs[0]);
	R = (double*)mxGetPr(prhs[0]);
	
	irodrigues(w,R);
	/* negera, ty R var col. major */
	w[0] = -w[0]; w[1] = -w[1]; w[2] = -w[2];
}


