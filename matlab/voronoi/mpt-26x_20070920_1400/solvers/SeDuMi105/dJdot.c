/*
 dJdot.c
 Jos F. Sturm, 2001.

   dJdotX = dJdot(d,X,K)

Given N x m matrix X, creates length(K.q) x m matrix dJdotX, having entries
d[i]'* J * xj[i] for each Lorentz block.
*/

#include <string.h>
#include "mex.h"
#include "blksdp.h"

#define DJDOTX_OUT plhs[0]
#define NPAROUT 1

#define D_IN prhs[0]
#define X_IN prhs[1]
#define K_IN prhs[2]
#define NPARIN 3

void dJdotxj(double *ypr, const double *d, const double *xpr,
            const int *lorNL, const int lorN)
{
  int k,nk;
  double yk;
  for(k = 0; k < lorN; k++){
    nk = lorNL[k];
    yk = d[0]*xpr[0];
    ypr[k] = yk - realdot(++d,++xpr, --nk);
    d += nk; xpr += nk;
  }
}


/* ============================================================
   MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
                 const int nrhs, const mxArray *prhs[])
{
  int lend, i, j, m;
  const double *d;
  double *dJdotxpr, *xpr;
  int *iwork;
  coneK cK;

/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  if(nrhs < NPARIN)
    mexErrMsgTxt("dJdot requires more input arguments.");
  if (nlhs > NPAROUT)
    mexErrMsgTxt("dJdot generates less output arguments.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Get statistics of cone K structure
   ------------------------------------------------------------ */
  lend = cK.lpN +  cK.qDim;
/* ------------------------------------------------------------
   Allocate working array iwork(cK.lorN+1).
   ------------------------------------------------------------ */
  iwork = (int *) mxCalloc(cK.lorN + 1, sizeof(int));
/* ------------------------------------------------------------
   Get scale data: d.
   ------------------------------------------------------------ */
  d = mxGetPr(D_IN);
  if(mxGetM(D_IN) * mxGetN(D_IN) != cK.qDim)
    if(mxGetM(D_IN) * mxGetN(D_IN) < lend)
      mexErrMsgTxt("d size mismatch");
    else
      d += cK.lpN;  /* point to lorentz */
/* ------------------------------------------------------------
   Get X
   ------------------------------------------------------------ */
  m = mxGetN(X_IN);
  if(mxGetM(X_IN) < lend)
    mexErrMsgTxt("X size mismatch");
  xpr = mxGetPr(X_IN);
  if(mxIsSparse(X_IN))
    mexErrMsgTxt("sparse x not supported");
  xpr += cK.lpN;                 /* point to Lorentz */
/* ------------------------------------------------------------
   If X is not sparse, then DJDOTX is full |K.q| x m. Compute for
   each column j = 1:m:
   ------------------------------------------------------------ */
  DJDOTX_OUT = mxCreateDoubleMatrix(cK.lorN, m, mxREAL);
  dJdotxpr = mxGetPr(DJDOTX_OUT);
/* ------------------------------------------------------------
   Let iwork = cK.lorNL in ints
   ------------------------------------------------------------ */
  for(i = 0; i < cK.lorN; i++)
    iwork[i] = cK.lorNL[i];
/* ------------------------------------------------------------
   Compute d[k]'*J*x[k,i] for all Lorentz blocks k.
   ------------------------------------------------------------ */
  lend = mxGetM(X_IN);           /* length of columns in X */
  for(i = 0; i < m; i++){
    dJdotxj(dJdotxpr, d, xpr, iwork, cK.lorN);
    dJdotxpr += cK.lorN;
    xpr += lend;                /* to next column */
  }
/* ------------------------------------------------------------
   Release common working arrays.
   ------------------------------------------------------------ */
  mxFree(iwork);
}
