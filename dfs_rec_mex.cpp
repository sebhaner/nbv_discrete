#include "mex.h"
#include <vector>

#define MAXDEPTH 2000
using namespace std;

mwIndex *ir, *jc;
mxArray *mx_nn;
int source, dest, nodes;
vector<int> pred;

bool search(const int node, const int depth)
{
	if(node>=nodes) mexErrMsgTxt("node out of range");
	if(node==dest) return true;
	int num_neigh = jc[node+1]-jc[node];
	if(num_neigh==0) return false;
	if(depth>MAXDEPTH){ mexPrintf("Warning, maximum recursion limit reached\n"); return false; }
	
	*mxGetPr(mx_nn) = num_neigh;
	mxArray *mx_rperm;
	mexCallMATLAB(1,&mx_rperm,1,&mx_nn,"randperm");
	double *rperm = mxGetPr(mx_rperm);
	
	for(int i=0; i<num_neigh; i++)
	{
		int neigh = ir[jc[node]+(int)rperm[i]-1];
		if(pred[neigh]==-1)
		{
			pred[neigh] = node;
			if(search(neigh,depth+1)) return true;
			pred[neigh] = -2; // mark as visited, do not search again (takes too long)
		}
	}
	mxDestroyArray(mx_rperm);
	return false;
}

// path = dfs_rec(Gt,source,dest)
void mexFunction(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
{
	
	if(nrhs!=3) mexErrMsgTxt("Requires 3 input arguments");
	source = (int)*mxGetPr(prhs[1]) - 1;
	dest = (int)*mxGetPr(prhs[2]) - 1;
	nodes = mxGetM(prhs[0]);
	if(mxGetN(prhs[0])!=nodes || !mxIsSparse(prhs[0])) mexErrMsgTxt("G must be square and sparse");
	
	if(source<0 || source>=nodes) mexErrMsgTxt("source out of range");
	if(dest<0 || dest>=nodes) mexErrMsgTxt("dest out of range");
	
	ir = mxGetIr(prhs[0]); // row indices
	jc = mxGetJc(prhs[0]); // column pointer
	
	mx_nn = mxCreateDoubleScalar(1);
// 	mexPrintf("%f\n",*mxGetPr(mx_nn));
	
	pred.resize(nodes);
	for(int i=0; i<nodes; i++) pred[i] = -1;
// 	pred[source] = -1;
	
	if(!search(source,0)) mexErrMsgTxt("Did not reach destination (search returned false)");
	if(pred[dest]<0) mexErrMsgTxt("Did not reach destination (pred[dest]<0)");
	
// 	for(int i=0; i<pred.size(); i++) mexPrintf("%d ",pred[i]);
// 	mexPrintf("\n");
	
	vector<int> path(nodes);
	path[0] = dest;
	int i;
	for(i=1; i<nodes; i++)
	{
		path[i] = pred[path[i-1]];
		if(path[i]==source) break;
	}
	
	
	plhs[0] = mxCreateDoubleMatrix(1,i+1,mxREAL);
	double *out = mxGetPr(plhs[0]);
	for(int j=0; j<=i; j++)
		out[j] = path[i-j]+1;
}





