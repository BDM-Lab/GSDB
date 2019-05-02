/*****************************************************************************
*   [HouseA, d] = mexhouse(A);
*   mex -O -largeArrayDims mexhouse.c
*   A \in R^{m*n} = QR,  QR Factorization via Householder vector
*   Assume m>= n,
*   each vector vk is stored in A(k:m,k), diag(R) is stored in the vector d
*   and triu(R,1) = triu(HouseA,1);
*   Algorithm 10.1, Numerical Linear Algegra
*   v4, July 15, 2009
* The performance of mexhousetest.c is nearly the same as mexhouse.c,
*  but this one is better.
*****************************************************************************/

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

/********************************
* realdotde: x dense matrix,  
*            y dense vector
*********************************/
double realdotde(const double *x, const int idx, 
                 const double *y, const int n)

{  int i;
   double r; 

   r=0.0;
   for (i=0; i<n-3; i++) {             /* LEVEL 4 */
      r += x[i+idx] * y[i]; i++;
      r += x[i+idx] * y[i]; i++;
      r += x[i+idx] * y[i]; i++;
      r += x[i+idx] * y[i]; }
   if (i<n-1) {                        /* LEVEL 2 */
      r += x[i+idx] * y[i]; i++;
      r += x[i+idx] * y[i]; i++; }
   if (i<n) {                          /* LEVEL 1 */
      r += x[i+idx] * y[i]; }  
   
   return r; 
}

double twonorm(const double *x, const int n)
{    int i;
      double sum, tmp;
      
      sum = 0.0;
      for (i=0; i<n; i++)
      {  
          tmp = x[i];
          if (tmp != 0)
          { 
              sum += tmp * tmp;
          } 
      }
       sum = sqrt(sum);
       
       return sum;
}


/**********************************************************/
void mexFunction(int nlhs,   mxArray  *plhs[], 
                 int nrhs,   const mxArray  *prhs[] )

{    double  *A, *d;
      
      double     *v;
      mxArray  *array_v;
      
      int i, j, k;
      int m, n, mk, m1k, nk, jm, j1m, j2m;
      double  alpha, tmp, vnorm;
      
      if (nrhs !=1) {
            mexErrMsgTxt(" mexAV2: requires only 1 input argument"); }
      if (nlhs >2) { 
            mexErrMsgTxt("mexAV2: requires 2 output argument"); }
      
      m = mxGetM(prhs[0]);
      n =  mxGetN(prhs[0]);
      
      if (m<n){
          mexErrMsgTxt("mexAV2: assumes row_size >= column_size"); }
      
      array_v = mxCreateDoubleMatrix(m, 1, mxREAL);
      v = mxGetPr(array_v);
   
      plhs[0] = mxDuplicateArray(prhs[0]);
      A = mxGetPr(plhs[0]);
      plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
      d = mxGetPr(plhs[1]);
      
      for (k=0; k< n; k++)
      {
          mk = m-k;  m1k = (m+1)*k;
          vnorm = twonorm(&A[m1k], mk);
          tmp = A[m1k];
          if (tmp>=0)
          {
              A[m1k]  += vnorm;
          }
          else
          {
              A[m1k]  -=  vnorm;
          } 
          vnorm = twonorm(&A[m1k], mk);
          alpha = sqrt(2)/vnorm;
          for (j=0; j<mk; j++)
          {
                 v[j] = alpha * A[m1k+j];
                 if (j==0)
                 {
                     A[m1k] = tmp;
                 }
          }
          
           nk = n-k;
           for (j=0; j<nk; j++)
          {
              j1m = j*mk; j2m = j*m + m1k;
              tmp = realdotde(A, j2m, v, mk);
              for (i=0; i<mk; i++)
              { 
                  if ((i==0)&&(j==0))
                 { 
                     d[k] = A[i+j2m] - v[i]*tmp;
                 }
                 if (j==0)
                 {
                      A[i+j2m] = v[i];
                 }
                  else
                  {
                      A[i+j2m] -= v[i]*tmp;
                  }
              }
          }
      }
      
      mxDestroyArray(array_v);  
      return;
}


