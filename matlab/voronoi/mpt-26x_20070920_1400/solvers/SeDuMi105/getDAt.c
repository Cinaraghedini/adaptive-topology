/*
 getDAt.c
 Jos F. Sturm, 2000.

   DAtq = getDAt(At,Ablk, d,K)

Creates DAtq  length(K.q) x m with Ablk.q sparsity structure.

*/

#include <string.h>
#include "mex.h"
#include "blksdp.h"

#define DATQ_OUT plhs[0]
#define NPAROUT 1

#define AT_IN prhs[0]
#define ABLK_IN prhs[1]
#define D_IN prhs[2]
#define K_IN prhs[3]
#define NPARIN 4

/* ============================================================
   L O R E N T Z  -  PART
   ============================================================ */
/* ************************************************************
   PROCEDURE ddotxj - Compute y[k] = d_k'*xj_k for each nonzero
      block in xj. (LORENTZ only)
   INPUT
     yir - lists Lorentz blocks where y has nonzeros, yir[i]<=lorN.
       the length of this array is implied by x, and is returned.
     d - lpN + qDim scaling vector.
     xir, xpr - sparse matrix. We compute d[k]'*xj[k] for each lorentz block
       where the column xj has nonzeros.
     xjc0, xjc1 - xir[xjc0:xjc1-1] should contain the Lorentz-subscripts for
       the desired column xj. xjc0 should point to the start of the Lorentz
       part, i.e. after LP part and before PSD part.
     blkstart - length lorN+1 array. Lorentz block k has subscripts
       blkstart[k]:blkstart[k+1]-1.
   OUTPUT
     ypr - length ynnz (returned) double array. Gives d[k]'*xj[k] for
       each Lorentz block listed in yir. This should correspond to the
       nonzero block structure in xir[xjc0:xjc1-1].
   RETURNS - ynnz. Number of entries (<= lorN) written into ypr.
   ************************************************************ */
int ddotxj(double *ypr, const int *yir, const double *d,
           const int *xir, const double *xpr,
           const int xjc0, const int xjc1,
           const int *blkstart)
{
  int knz, nexti, inz, i;
  double yk;
/* ------------------------------------------------------------
   If x is all-0, then return 0.
   ------------------------------------------------------------ */
  if(xjc0 >= xjc1)
    return 0;
/* ------------------------------------------------------------
   INITIALIZE: knz points to nz's (Lorentz blocks) in y.
   The current block has subscripts smaller than nexti.
   Accumulate yk = ddotxj[k]. Init to 0.
   ------------------------------------------------------------ */
  knz = 0;
  nexti = blkstart[yir[knz] + 1];              /* yir lists Lorentz blocks */
  yk = 0.0;
/* ------------------------------------------------------------
   Browse through nonzeros in xj
   ------------------------------------------------------------ */
  for(inz = xjc0; inz < xjc1; inz++)
    if( (i = xir[inz]) < nexti)
      yk += d[i] * xpr[inz];
/* ------------------------------------------------------------
   If we finished the current nonzero Lorentz block, then write entry.
   ------------------------------------------------------------ */
    else{
      ypr[knz++] = yk;
      yk = d[i] * xpr[inz];
      nexti = blkstart[yir[knz] + 1];       /* yir lists Lorentz blocks */
    }
/* ------------------------------------------------------------
   Write last yk = ddotxj[k] entry into y, and return nnz(y).
   ------------------------------------------------------------ */
  ypr[knz++] = yk;
  return knz;
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
  const mxArray *MY_FIELD;
  int lenfull, lenud, i, j, m, inz, nblk;
  const int *blkstart;
  const double *d, *AblkjcPr, *blkstartPr;
  int *Ablkjc, *Ablkjc1, *blkstartOrg;
  coneK cK;
  jcir At, ddota;

/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  if(nrhs < NPARIN)
    mexErrMsgTxt("getDAt requires more input arguments.");
  if (nlhs > NPAROUT)
    mexErrMsgTxt("getDAt generates less output arguments.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Get statistics of cone K structure
   ------------------------------------------------------------ */
  lenud = cK.rDim + cK.hDim;    /* dim PSD part */
  lenfull = cK.lpN +  cK.qDim + lenud;
  nblk = 1 + cK.lorN + cK.sdpN;
/* ------------------------------------------------------------
   Allocate working array blkstart(nblk+1).
   ------------------------------------------------------------ */
  blkstartOrg = (int *) mxCalloc(nblk + 1, sizeof(int));
/* ------------------------------------------------------------
   Translate blkstart from Fortran-double to C-int
   ------------------------------------------------------------ */
  if( (MY_FIELD = mxGetField(K_IN,0,"blkstart"))==NULL)  /*K.blkstart*/
    mexErrMsgTxt("Missing K.blkstart.");
  if(mxGetM(MY_FIELD) * mxGetN(MY_FIELD) != nblk + 1)    /* check nblk+1 */
    mexErrMsgTxt("Size mismatch K.blkstart.");
  blkstartPr = mxGetPr(MY_FIELD);
  for(i = 0; i <= nblk; i++){                         /* to integers */
    j = blkstartPr[i];
    blkstartOrg[i] = --j;
  }
  blkstart = blkstartOrg;
/* ------------------------------------------------------------
   Get scale data: d.
   ------------------------------------------------------------ */
  if(mxGetM(D_IN) * mxGetN(D_IN) < lenfull - lenud)
    mexErrMsgTxt("d size mismatch");
  d = mxGetPr(D_IN);
/* ------------------------------------------------------------
   Get At
   ------------------------------------------------------------ */
  if(!mxIsSparse(AT_IN))
    mexErrMsgTxt("At must be sparse");
  m = mxGetN(AT_IN);
  if(mxGetM(AT_IN) != lenfull)
    mexErrMsgTxt("At size mismatch");
  At.jc = mxGetJc(AT_IN);
  At.ir = mxGetIr(AT_IN);
  At.pr = mxGetPr(AT_IN);
/* ------------------------------------------------------------
   Disassemble Ablk structure
   ------------------------------------------------------------ */
  if(!mxIsStruct(ABLK_IN))
    mexErrMsgTxt("Ablk should be a structure.");
  if( (MY_FIELD = mxGetField(ABLK_IN,0,"jc")) == NULL)      /* Ablk.jc */
    mexErrMsgTxt("Missing field Ablk.jc.");
  if(mxGetM(MY_FIELD) != m || mxGetN(MY_FIELD) != 2)
    mexErrMsgTxt("Ablk.jc size mismatch");
  AblkjcPr = mxGetPr(MY_FIELD);
  if( (MY_FIELD = mxGetField(ABLK_IN,0,"q")) == NULL)       /* Ablk.q */
    mexErrMsgTxt("Missing field Ablk.q.");
  if(mxGetM(MY_FIELD) != cK.lorN || mxGetN(MY_FIELD) != m)
    mexErrMsgTxt("Ablk.q size mismatch");
  if(!mxIsSparse(MY_FIELD))
    mexErrMsgTxt("Ablk.q must be sparse");
/* ------------------------------------------------------------
   Create output DAtq as copy of Ablk.q
   ------------------------------------------------------------ */
  DATQ_OUT = mxDuplicateArray(MY_FIELD);              /* init DAtq = Ablk.q */
  ddota.jc = mxGetJc(DATQ_OUT);
  ddota.ir = mxGetIr(DATQ_OUT);
  ddota.pr = mxGetPr(DATQ_OUT);
/* ------------------------------------------------------------
   Allocate working arrays:
   int Ablkjc(2*m)
   double fwork max(cK.rMaxn^2, 2*cK.hMaxn^2).
   ------------------------------------------------------------ */
  Ablkjc = (int *) mxCalloc(2*m, sizeof(int) );
/* ------------------------------------------------------------
   Convert double to int:
   ------------------------------------------------------------ */
  for(i = 0; i < 2 * m; i++){                /* Ablkjc */
    Ablkjc[i] = AblkjcPr[i];
  }
  Ablkjc1 = Ablkjc + m;
/* ------------------------------------------------------------
   LORENTZ:
   ------------------------------------------------------------ */
  ++blkstart;               /* point to Lorentz part */
  inz = 0;
  for(i = 0; i < m; i++){
    inz += ddotxj(ddota.pr + inz, ddota.ir + inz, d, At.ir, At.pr,
                  Ablkjc[i], Ablkjc1[i], blkstart);
    mxAssert(inz == ddota.jc[i+1],"");
  }
/* ------------------------------------------------------------
   Release working arrays.
   ------------------------------------------------------------ */
  mxFree(Ablkjc);
  mxFree(blkstartOrg);
}
