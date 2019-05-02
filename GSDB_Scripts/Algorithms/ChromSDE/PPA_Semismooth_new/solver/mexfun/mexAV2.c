/************************************************
*    AV2 = mexAV2(A, Q, options)
*    mex -O -largeArrayDims -lmwlapack -lmwblas mexAV2.c
*    options = 0 (default), A*V2;
*    options = 1, A*(V2t);
*    householder vectors v are stored in tril(Q);
*    v9, July 15, 2009
*************************************************/

#include <mex.h>
#include <math.h>
#include <matrix.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if !defined(MX_API_VER) || ( MX_API_VER < 0x07030000 )
typedef int mwIndex;
typedef int mwSize;
#endif

#if !defined(MAX)
#define  MAX(A, B)   ((A) > (B) ? (A) : (B))
#endif

void mexFunction(int nlhs,   mxArray  *plhs[], 
                 int nrhs,   const mxArray  *prhs[] )
{    
      double  *A, *Q, *lhspr, *rhspr;
      mxArray  *rhsvalue;
      
      double      *Av;
      mxArray    *array_Av;
      
      int         i, j, k, options;
      int         m, n, m1, n1, jm, j1m, j2m;
      double   tmp;
      
     mwSize   M, N, INCX, INCY;
     mwSize   LDA;
     double    ALPHA, BETA;
     char        *TRANS="N";

      if (nrhs <2) {
            mexErrMsgTxt(" mexAV2: must have at least 2 inputs"); }
      if (nlhs > 1) { 
            mexErrMsgTxt("mexAV2: requires 1 output argument"); }
      if (nrhs == 2) 
            { options = 0; } 
      else 
           { options = (int) *mxGetPr(prhs[2]); }
      
      m = mxGetM(prhs[0]);
      n =  mxGetN(prhs[0]);
      m1 = mxGetM(prhs[1]); 
      n1 = mxGetN(prhs[1]);
      Q = mxGetPr(prhs[1]);
      
     ALPHA = 1.0;
     BETA =0.0;
     INCX = 1;
     INCY = 1;

   if (options==0)
   {
      
      if (n != m1)
      {    
          mexErrMsgTxt("mexAV2:1ST and 2ND input not compatible.");
      }
      
      array_Av = mxCreateDoubleMatrix(m, 1, mxREAL);
      Av = mxGetPr(array_Av);
      
      rhsvalue = mxDuplicateArray(prhs[0]);
      A = mxGetPr(rhsvalue);
      
      plhs[0] = mxCreateDoubleMatrix(m, m1-n1, mxREAL);
      lhspr = mxGetPr(plhs[0]);
      M = m;
      LDA = M;


         for ( k=0; k<n1; k++)
            {    
                N = n-k;
                 
  dgemv_(TRANS, &M, &N, &ALPHA, &A[m*k], &LDA, &Q[(m1+1)*k], &INCX, &BETA, Av, &INCY);
 
                 for( j=0; j< N; j++)
                 { 
                     j1m = j*m;  j2m = j*m + m*k;
                     tmp = Q[(m1+1)*k+ j];
                     for( i=0; i<m; i++)
                     {  
                         A[i+j2m] -= Av[i]*tmp;
                     }
                 }
            }
        /***************************************/
         for ( j=0; j<m*(m1-n1); j++)
         {  
             lhspr[j] = A[m*n1+j];
         }
      
         mxDestroyArray(rhsvalue);
   }
   else
   {   
       if (n != (m1-n1) )
          mexErrMsgTxt("mexAV2:1ST and 2ND input not compatible.");
       
       plhs[0] = mxCreateDoubleMatrix(m, m1, mxREAL);
       A = mxGetPr(plhs[0]);
       rhspr = mxGetPr(prhs[0]);
       for ( j=0; j< m*(m1-n1); j++)
       {    A[m*n1+j] =  rhspr[j];    }
     
      array_Av = mxCreateDoubleMatrix(m, 1, mxREAL);  
      Av = mxGetPr(array_Av);
      M = m;
      LDA = M;

      for ( k=(n1-1); k>=0; k--)
            {
               N = m1-k;
                 
dgemv_(TRANS, &M, &N, &ALPHA, &A[m*k], &LDA, &Q[(m1+1)*k], &INCX, &BETA, Av, &INCY);
               
               for( j=0; j< N; j++)
                 { 
                     j1m = j*m;  j2m = j*m + m*k;
                     tmp = Q[(m1+1)*k+ j];
                     for( i=0; i<m; i++)
                     {  
                         A[i+j2m] -= Av[i] * tmp;
                     }
                 }
            }
     }
      
      mxDestroyArray(array_Av);
      
      return;
}

