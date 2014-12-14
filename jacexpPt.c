#include <math.h>
#include "mex.h"

#define THRESH (0.1e-6)

void JacPt(double A[3][2], double w1, double w2, double w3,
		double a1, double a2, double a3, double X, double Y, double Z){
	
	double t10, t14, t15, t17, t18, t19, t20, t21, t22, t25, t27, t28,
			t3, t33, t35, t36, t4, t40, t42, t43, t45, t48, t49, t5,
			t50, t56, t6, t62, t64, t67, t7, t8, t9;
	
	t3 = w1*w1;
	t4 = w2*w2;
	t5 = w3*w3;
	t6 = t3 + t4 + t5;
	t7 = sqrt(t6);
	t8 = cos(t7);
	t9 = 0.1e1 - t8;
	t10 = 0.1e1 / t6;
	t14 = 0.1e1 + t9 * (t3 * t10 - 0.1e1);
	t15 = sin(t7);
	t17 = 0.1e1 / t7;
	t18 = t15 * w2 * t17;
	t19 = t9 * w1;
	t20 = t10 * w3;
	t21 = t19 * t20;
	t22 = -t18 + t21;
	t25 = t15 * w1 * t17;
	t27 = t9 * w2 * t20;
	t28 = t25 + t27;
	t33 = 0.1e1 + t9 * (t5 * t10 - 0.1e1);
	t35 = t22 * X + t28 * Y + t33 * Z + a3;
	t36 = 0.1e1 / t35;
	t40 = t15 * w3 * t17;
	t42 = t19 * t10 * w2;
	t43 = -t40 + t42;
	t45 = t18 + t21;
	t48 = t35 * t35;
	t49 = 0.1e1 / t48;
	t50 = (t14 * X + t43 * Y + t45 * Z + a1) * t49;
	A[0][0] = t14 * t36 - t50 * t22;
	A[1][0] = t43 * t36 - t50 * t28;
	A[2][0] = t45 * t36 - t50 * t33;
	t56 = t40 + t42;
	t62 = 0.1e1 + t9 * (t4 * t10 - 0.1e1);
	t64 = -t25 + t27;
	t67 = (t56 * X + t62 * Y + t64 * Z + a2) * t49;
	A[0][1] = t56 * t36 - t67 * t22;
	A[1][1] = t62 * t36 - t67 * t28;
	A[2][1] = t64 * t36 - t67 * t33;
}



void JacPtSmallW(double A[3][2], double w1, double w2, double w3,
		double a1, double a2, double a3, double X, double Y, double Z){
	
	double t10, t11, t12, t14, t20, t5, t6;
	
	t5 = -w2 * X + w1 * Y + Z + a3;
	t6 = 0.1e1 / t5;
	t10 = t5 * t5;
	t11 = 0.1e1 / t10;
	t12 = (X - w3 * Y + w2 * Z + a1) * t11;
	A[0][0] = t6 + t12 * w2;
	t14 = w3 * t6;
	A[1][0] = -t14 - t12 * w1;
	A[2][0] = w2 * t6 - t12;
	t20 = (w3 * X + Y - w1 * Z + a2) * t11;
	A[0][1] = t14 + t20 * w2;
	A[1][1] = t6 - t20 * w1;
	A[2][1] = -w1 * t6 - t20;
}

double *out1, *w, *a, *X;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	if(nrhs!=3) mexErrMsgTxt("Needs 3 input arguments, each of length 3");
	if(nlhs>1) mexErrMsgTxt("Too many output variables");
	w = mxGetPr(prhs[0]);
	a = mxGetPr(prhs[1]);
	X = mxGetPr(prhs[2]);
	
	plhs[0] = mxCreateDoubleMatrix(2, 3, 0);
	out1 = mxGetPr(plhs[0]);
	
	if(w[0]*w[0]+w[1]*w[1]+w[2]*w[2] > THRESH*THRESH)
		JacPt(out1, w[0], w[1], w[2], a[0], a[1], a[2], X[0], X[1], X[2]);
	else JacPtSmallW(out1, w[0], w[1], w[2], a[0], a[1], a[2], X[0], X[1], X[2]);
	
	
}