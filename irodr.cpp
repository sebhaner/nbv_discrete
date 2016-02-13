#include "mex.h"
#include "rodrigues.cpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	if(nrhs!=1) mexErrMsgTxt("Needs 1 input argument");
	if(!mxIsDouble(prhs[0]) && !mxIsSingle(prhs[0])) mexErrMsgTxt("Input must be of class DOUBLE or SINGLE");
	if(mxGetNumberOfElements(prhs[0])!=9) mexErrMsgTxt("Input must have 9 elements");
	
	if(mxIsDouble(prhs[0]))
	{
		
		plhs[0] = mxCreateDoubleMatrix(3,1,mxREAL);
		double* w = (double*)mxGetPr(plhs[0]);
		double* R = (double*)mxGetPr(prhs[0]);
	
		irodrigues(w,R);
		/* negera, ty R var col. major */
		w[0] = -w[0]; w[1] = -w[1]; w[2] = -w[2];
		
	} else {
		
		float *Rf = (float*)mxGetPr(prhs[0]);
		double R[9], w[3];
		R[0] = Rf[0]; R[1] = Rf[3]; R[2] = Rf[6];
		R[3] = Rf[1]; R[4] = Rf[4]; R[5] = Rf[7];
		R[6] = Rf[2]; R[7] = Rf[5]; R[8] = Rf[8];
		irodrigues(w,R);
		
		plhs[0] = mxCreateNumericMatrix(3,1,mxSINGLE_CLASS,mxREAL);
		float* wf = (float*)mxGetPr(plhs[0]);
		wf[0] = w[0]; wf[1] = w[1]; wf[2] = w[2];
		
	}
}


