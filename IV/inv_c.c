#include "mex.h"

void inv(mxDouble *A, mxDouble *outM, mwSize n)
{
  mwSize i, j, k, offset;
  mxDouble tmp;
  
  for(i = 0; i < n; i++)      /* outM - unit matrix */ 
    for(j = 0; j < n; j++)
      if(i == j)
        outM[j+i*n] = 1.0;
  
  for(k = 0; k < n; k++){
    tmp = A[k+k*n];      
    for(j = 0; j < n; j++){   /* M[k][j] = M[k][j]/M[k][k] */
      offset = j + k*n;
      A[offset] /= tmp;
      outM[offset] /= tmp;
    }
    for(i = 0; i < n; i++){   
      tmp = A[k+i*n];
      for(j = 0; j < n; j++){ /* M[i][j] = M[i][j] - M[k][j]*M[i][k] */
        if(i == k)
          break;
        offset = j + i*n;
        A[offset]    -= A[j+k*n]*tmp;
        outM[offset] -= outM[j+k*n]*tmp;
      }
    }
  }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  mxDouble *A, *outM;
  mwSize n;
  
  if(nlhs > 1)
    mexErrMsgIdAndTxt("IV:inv_c:nlhs", "Expected single output.\n");
  if(nrhs != 1)
    mexErrMsgIdAndTxt("IV:inv_c:nrhs", "Expected single input.\n");
  if(mxGetN(prhs[0]) != mxGetM(prhs[0]))
    mexErrMsgIdAndTxt("IV:inv_c:prhs", "Matrix is not square.\n");
  /* Assume non-zero determinant */
  
  n = (mwSize)mxGetN(prhs[0]);
  A = mxGetDoubles(prhs[0]);
  
  plhs[0] = mxCreateDoubleMatrix(n,n,mxREAL);
  outM = mxGetDoubles(plhs[0]);
  inv(A,outM,n);
}