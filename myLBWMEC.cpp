#include "mex.h" 
#include "mytest.h"

// function type
// [Y,V,F] = myLBWMEC(X,P)
// size of each matrix 
// Y: M*2, V:M*1, F:M*1, X: N*2, P: M*2

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] ) 
// the first two are input, the last three are output
{ 
    // check parameter number
    if(nrhs!=2){ 
        mexErrMsgTxt("One input required."); 
    }else if(nlhs>3){ 
        mexErrMsgTxt("Too many output arguments"); 
    } 
    // get the size of input matrix    
    int Xn = mxGetN(prhs[0]); 
    int Pn = mxGetN(prhs[1]);
    
    // construct output matrix
    plhs[0] = mxCreateDoubleMatrix(2,Pn,mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(1,Pn,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,Pn,mxREAL);
    
    // get pointer    
    double *X = mxGetPr(prhs[0]); 
    double *P = mxGetPr(prhs[1]);
    double *Y = mxGetPr(plhs[0]); 
    double *V = mxGetPr(plhs[1]);
    double *F = mxGetPr(plhs[2]);
    
    // call C++ function
    myTest test(X,P,Xn,Pn);
    test.buildTreeStructure();
    test.findsolution(Y,V,F);
} 