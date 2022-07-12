#include <math.h>

#include "mex.h"
#include "matrix.h"

mxComplexDouble sqrtComplex(mxComplexDouble D)
{
  mxComplexDouble Dsqrt;
  double Dabs, Darg;
  
  Dabs = sqrt(pow(D.real, 2.0)+pow(D.imag, 2.0));
  Darg = atan2(D.imag, D.real);
  
  Dsqrt.real = sqrt(Dabs)*cos(Darg/2);
  Dsqrt.imag = sqrt(Dabs)*sin(Darg/2);
  
  return Dsqrt;
}

mxComplexDouble getRoot(mxComplexDouble A, mxComplexDouble B, 
                        mxComplexDouble Dsqrt)
{
  mxComplexDouble X;
  double Aabs2;
  
  Aabs2 = pow(A.real,2.0)+pow(A.imag,2.0);

  X.real =  (A.real/Aabs2)*(-B.real + Dsqrt.real)/2.0 + 
            (A.imag/Aabs2)*(-B.imag + Dsqrt.imag)/2.0;
  
  X.imag =  (A.real/Aabs2)*(-B.imag + Dsqrt.imag)/2.0 -
            (A.imag/Aabs2)*(-B.real + Dsqrt.real)/2.0;
  return X;
}

mxComplexDouble getDisc(mxComplexDouble A, mxComplexDouble B,
                        mxComplexDouble C)
{
  mxComplexDouble D;

  D.real = pow(B.real,2.0) - pow(B.imag,2.0) + 
            4*(A.imag*C.imag - A.real*C.real);
  D.imag = 2*B.real*B.imag - 
            4*(A.real*C.imag + A.imag*C.real);
  return D;
}

/* The computational routine */
void quadsolve( mxComplexDouble *A, mxComplexDouble *B, mxComplexDouble *C, 
                mxArray *X1, mxArray *X2, mxArray *D, 
                mwSize rows, mwSize cols)
{
  mwSize i, j, offset;
  mxComplexDouble * X1c = mxGetComplexDoubles(X1);
  mxComplexDouble * X2c = mxGetComplexDoubles(X2);
  mxComplexDouble * Dc  = mxGetComplexDoubles(D);
  mxComplexDouble Dsqrt;
  
  for(i = 0; i < rows; i++)
    for(j = 0; j < cols; j++){
      offset = j + rows*i;
      
      Dc[offset] = getDisc(A[offset],B[offset],C[offset]);
      Dsqrt = sqrtComplex(Dc[offset]);
      
      X1c[offset] = getRoot(A[offset],B[offset],Dsqrt);
      
      Dsqrt.real  = -Dsqrt.real;
      Dsqrt.imag  = -Dsqrt.imag;
      
      X2c[offset] = getRoot(A[offset],B[offset],Dsqrt);
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mxComplexDouble *coeff_A;            
    mxComplexDouble *coeff_B;
    mxComplexDouble *coeff_C;
    mwSize rows;                     
    mwSize cols;
    mxArray *tmp;         /* pass as D when it is not required */

    /* check for proper number of arguments */
    if(nrhs != 3)
      mexErrMsgIdAndTxt("IV:quadsolve:nrhs",
                        "Three shall be the number of arguments thou shalt input.\n");
    if(nlhs < 2)
      mexErrMsgIdAndTxt("IV:quadsolve:nlhs",
                        "At least 2 outputs required.\n");
    if(nlhs > 3)
      mexErrMsgIdAndTxt("IV:quadsolve:nlhs",
                        "No more than 3 outputs required.\n");
    /* Check that all inputs are complex*/
    if(!mxIsComplex(prhs[0]) || !mxIsComplex(prhs[1]) || !mxIsComplex(prhs[2]))
      mexErrMsgIdAndTxt("IV:quadsolve:inputsNotComplex",
                        "Inputs must be complex.\n");
    
    rows = (mwSize)mxGetM(prhs[0]);
    cols = (mwSize)mxGetN(prhs[0]);
    
    coeff_A = mxGetComplexDoubles(prhs[0]);
    coeff_B = mxGetComplexDoubles(prhs[1]);
    coeff_C = mxGetComplexDoubles(prhs[2]);
    
    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxCOMPLEX);
    plhs[1] = mxCreateDoubleMatrix(rows, cols, mxCOMPLEX);
    
    if(nlhs == 3){ /* D is required argument */
      plhs[2] = mxCreateDoubleMatrix(rows, cols, mxCOMPLEX);
      quadsolve(coeff_A,coeff_B,coeff_C,plhs[0],plhs[1],plhs[2],rows,cols);
    } else {
      tmp = mxCreateDoubleMatrix(rows, cols, mxCOMPLEX);
      quadsolve(coeff_A,coeff_B,coeff_C,plhs[0],plhs[1],tmp,rows,cols);
    }
}
