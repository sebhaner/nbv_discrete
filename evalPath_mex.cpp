#include "mex.h"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

#define MAX(a,b)((a)>(b) ? (a) : (b))
#define MIN(a,b)((a)<(b) ? (a) : (b))

using namespace Eigen;

bool binsearch(const mwIndex *ir, int low, int high, mwIndex row, mwIndex &res){
	mwIndex m;
	while(low<high){
		m = (low + high)/2;
		if(ir[m] > row) high = m;
		else if(ir[m] < row) low = m + 1;
		else{ low = m; break; }
	}
	if(ir[low]==row){ res = low; return true; }
	else return false; // not found
}

// val = evalPath(path,G,info,lambda,type,[penalty_weight])
void mexFunction(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
{
	
	if(nrhs<5 || nrhs>6) mexErrMsgTxt("Requires 5 or 6 input arguments");
	const int plen = mxGetNumberOfElements(prhs[0]);
	const double *path = mxGetPr(prhs[0]);
	const int gm = mxGetM(prhs[1]);
	if(mxGetN(prhs[1])!=gm || !mxIsSparse(prhs[1])) mexErrMsgTxt("G must be square and sparse");
	if(!mxIsCell(prhs[2])) mexErrMsgTxt("info must be a cell array");
	if(mxGetNumberOfElements(prhs[2])!=gm+1) mexErrMsgTxt("info/G size mismatch");
	
	const int ndim = mxGetNumberOfDimensions(mxGetCell(prhs[2],0));
	const mwSize *dims = mxGetDimensions(mxGetCell(prhs[2],0));
	int nop;
	if(ndim!=3) nop = 1;
	else nop = dims[2];
	if(nop==0) mexErrMsgTxt("info must contain at least one point");
	
	if(mxIsEmpty(prhs[3])||!mxIsDouble(prhs[3])) mexErrMsgTxt("param should be a double scalar");
	const double lambda = *mxGetPr(prhs[3]);
	
	const int type = (int)*mxGetPr(prhs[4]);
	if(type<1 || type>4) mexErrMsgTxt("type must be 1, 2, 3 or 4");
	
	double penalty_weight = 1e6;
	if(nrhs==6){
		if(mxIsEmpty(prhs[5])||!mxIsDouble(prhs[5]))
			mexErrMsgTxt("penalty_weight should be a double scalar");
		else penalty_weight = *mxGetPr(prhs[5]);
	}
	
	const mwIndex *ir = mxGetIr(prhs[1]); // row indices
	const mwIndex *jc = mxGetJc(prhs[1]); // column pointer
	const double *gw = mxGetPr(prhs[1]);
	
	double len = 0;
	mwIndex val_ind;
	for(int i=0; i<plen-1; i++)
	{
		// find G(path(i),path(i+1))
		int col = path[i+1]-1;
		int row = path[i]-1;
		if(row>=gm || col>=gm) mexErrMsgTxt("Path index out of range");
		if(!binsearch(ir,jc[col],jc[col+1]-1,row,val_ind)) mexErrMsgTxt("Invalid path (or other error)");
		len += gw[val_ind];
	}
// 	mexPrintf("len: %f, nop %u, ndim: %u\n",len,nop,ndim);
	
	double geo = 0;
	for(int p=0; p<nop; p++)
	{
		Map<Matrix3d> I0(mxGetPr(mxGetCell(prhs[2],0))+9*p,3,3);
		Matrix3d Isum;
		Isum = I0;
		for(int k=0; k<plen; k++)
		{
			Map<Matrix3d> I(mxGetPr(mxGetCell(prhs[2],path[k]))+9*p,3,3);
			Isum += I; 
		}
		if(type<=2)
			geo += Isum.inverse().trace();
		else
		{
			double maxeig = 1.0/Isum.selfadjointView<Lower>().eigenvalues().minCoeff();
			if(maxeig>geo) geo = maxeig;
		}
	}
	
	if(type==1 || type==3)
		plhs[0] = mxCreateDoubleScalar(len+geo/lambda);
	else
	{
		double penalty = penalty_weight*MAX(0,len-lambda)*MAX(0,len-lambda);
		plhs[0] = mxCreateDoubleScalar(geo+penalty);
	}
	if(nlhs>1) plhs[1] = mxCreateDoubleScalar(len);
	if(nlhs>2) plhs[2] = mxCreateDoubleScalar(geo);
}





