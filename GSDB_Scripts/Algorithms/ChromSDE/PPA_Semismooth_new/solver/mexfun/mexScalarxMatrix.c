/************************************************************************
*  mexScalarxMatrix(X,alpha,options)         
*  X = alpha*X; 
************************************************************************/

#include "mex.h"
#include <math.h>
#include <matrix.h>

#if !defined(MX_API_VER) || ( MX_API_VER < 0x07030000 )
typedef int mwIndex;
typedef int mwSize;
#endif

/********************************************************************
  PROCEDURE mexFunction - Entry for Matlab
*********************************************************************/
void mexFunction(const int nlhs, mxArray *plhs[],
                 const int nrhs, const mxArray *prhs[])
{
  double   *X;
  mwIndex  *irX, *jcX; 
  int       m, n, isspX, j, jm, k, kstart, kend; 
  double    alpha;

  if(nrhs < 2)
    mexErrMsgTxt("mexschurfun: requires at least 2 input arguments.");
  if(nlhs > 0)
    mexErrMsgTxt("mexschurfun: requires no output argument.");

  X = mxGetPr(prhs[0]);
  isspX = mxIsSparse(prhs[0]); 
  if (isspX) {
     irX = mxGetIr(prhs[0]);
     jcX = mxGetJc(prhs[0]);
  }
  m = mxGetM(prhs[0]); 
  n = mxGetN(prhs[0]); 
  alpha = *mxGetPr(prhs[1]); 
  /********************************************************/
     if (isspX) { 
        for (j=0; j<n; j++) {
	   kstart = jcX[j]; kend = jcX[j+1]; 
           for (k=kstart; k<kend; k++) { 
              X[k] *= alpha; }
	}
     } else {
        for (j=0; j<n; j++) { 
           jm = j*m;
           for (k=0; k<m; k++) { X[k+jm] *= alpha; }
	}
     }   
return;
}
/************************************************************************/
