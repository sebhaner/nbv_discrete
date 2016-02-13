#include <math.h>
#define PI (3.141592653589793)
#define TH (1e-6)
#define MAX(a,b)((a)>(b) ? (a):(b))

void rodrigues(double *R, double w1, double w2, double w3){
	
	double theta = sqrt(w1*w1+w2*w2+w3*w3);
	if(theta<TH){
		R[0] = 1;   R[1] = -w3; R[2] = w2;
		R[3] = w3;  R[4] = 1;   R[5] = -w1;
		R[6] = -w2; R[7] = w1;  R[8] = 1;
	}
	else{
		w1 /= theta; w2 /= theta; w3 /= theta;
		double cm1 = cos(theta) - 1, s = sin(theta);
		double w1w2c = w1*w2*cm1, w1w3c = w1*w3*cm1, w2w3c = w2*w3*cm1;
		
		R[0] = 1-(w1*w1-1)*cm1;  R[1] = -w3*s-w1w2c;	 R[2] = w2*s-w1w3c;
		R[3] = w3*s-w1w2c;		 R[4] = 1-(w2*w2-1)*cm1; R[5] = -w1*s-w2w3c;
		R[6] = -w2*s-w1w3c;		 R[7] = w1*s-w2w3c;		 R[8] = 1-(w3*w3-1)*cm1;
	}
}

void irodrigues(double *w, double *R){
	
	double theta = (R[0]+R[4]+R[8]-1)*0.5;
	if(theta<-1) theta = -1; // gardera mot avrundningsfel i R
	if(theta>1) theta = 1;
	theta = acos(theta);
	
	if(theta>(PI-TH)){ // singulärt nära pi, använd ww' = (R-I)/2
// 		w[0] = sqrt((R[0]+1)*0.5)*theta; // ger abs(w)
// 		w[1] = sqrt((R[4]+1)*0.5)*theta;
// 		w[2] = sqrt((R[8]+1)*0.5)*theta;
		w[0] = (R[0]+1)*0.5;
		w[1] = (R[4]+1)*0.5;
		w[2] = (R[8]+1)*0.5;
		w[0] = sqrt(MAX(w[0],0))*theta; // i sällsynta fall kan det bli negativt
		w[1] = sqrt(MAX(w[1],0))*theta;
		w[2] = sqrt(MAX(w[2],0))*theta;
		
		if(w[0]>=w[1] && w[0]>=w[2]){ // sätt största elt. +, härled övriga
			w[1] = R[1]>=0 ? w[1] : -w[1];
			w[2] = R[2]>=0 ? w[2] : -w[2];
		}else if(w[1]>=w[0] && w[1]>=w[2]){
			w[0] = R[1]>=0 ? w[0] : -w[0];
			w[2] = R[5]>=0 ? w[2] : -w[2];
		}else{
			w[0] = R[2]>=0 ? w[0] : -w[0];
			w[1] = R[5]>=0 ? w[1] : -w[1];
		}
	}
	else{
		w[0] = R[7]-R[5]; w[1] = R[2]-R[6]; w[2] = R[3]-R[1];
		if(theta>TH){
			theta /= 2*sin(theta);
			w[0] *= theta; w[1] *= theta; w[2] *= theta;
		}
		else{
			w[0] *= 0.5; w[1] *= 0.5; w[2] *= 0.5;
		}
	}
}

