#include "mex.h"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

using namespace Eigen;

// val = evalGeoCost(path,info,type)
void mexFunction(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
{
	
	if(nrhs!=3) mexErrMsgTxt("Requires 3 input arguments");
	const int plen = mxGetNumberOfElements(prhs[0]);
	const double *path = mxGetPr(prhs[0]);
	if(!mxIsCell(prhs[1])) mexErrMsgTxt("info must be a cell array");
	
	const int ndim = mxGetNumberOfDimensions(mxGetCell(prhs[1],0));
	const mwSize *dims = mxGetDimensions(mxGetCell(prhs[1],0));
	int nop;
	if(ndim!=3) nop = 1;
	else nop = dims[2];
	if(nop==0) mexErrMsgTxt("info must contain at least one point");
	
	const int nodes = mxGetNumberOfElements(prhs[1])-1;
	
	const int type = (int)*mxGetPr(prhs[2]);
	if(type<1 || type>4) mexErrMsgTxt("type must be 1, 2, 3 or 4");
	
	
	double geo = 0;
	for(int p=0; p<nop; p++)
	{
		Map<Matrix3d> I0(mxGetPr(mxGetCell(prhs[1],0))+9*p,3,3);
		Matrix3d Isum;
		Isum = I0;
		for(int k=0; k<plen; k++)
		{
			if(path[k]>nodes) mexErrMsgTxt("Path index/info size mismatch");
			Map<Matrix3d> I(mxGetPr(mxGetCell(prhs[1],path[k]))+9*p,3,3);
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
	

	plhs[0] = mxCreateDoubleScalar(geo);

}





