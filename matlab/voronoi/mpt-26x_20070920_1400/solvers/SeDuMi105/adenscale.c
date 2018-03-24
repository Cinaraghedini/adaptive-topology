/*
 smult = adenscale(dense, d, qblkstart)

    This file is part of SeDuMi 1.05 TEST
    Copyright (C) 2001 Jos F. Sturm
    Dept. Quantitative Economics, Maastricht University, the Netherlands.
    Affiliations up to SeDuMi 1.02 (AUG1998):
      CRL, McMaster University, Canada.
      Supported by the Netherlands Organization for Scientific Research (NWO).

  
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
  
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
  
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/
#include <math.h>
#include <string.h>
#include "mex.h"
#include "blksdp.h"

#define SMULT_OUT plhs[0]
#define NPAROUT 1

#define DENSE_IN prhs[0]
#define D_IN prhs[1]
#define BLKSTART_IN prhs[2]
#define NPARIN 3

/* ************************************************************
   PROCEDURE adenscale - Computed Lorentz norm-bound part of smult such that
     AP(d)A' = ADA + Ad*diag(smult)*Ad'. Thus, smult(j) = det(k) if j
     belongs to Lorentz block k.
   INPUT
     detd     - length |K.q| vector of det(dk)_k for Lorentz; use its sqrt.
     dencols  - length nden array
     q        - length nq array
     blkend   - length nq array, with 1-beyond subscript for each Lorentz
        block listed in q.
     nq       - number of (removed) dense Lorentz blocks
     nden     - number of dense columns, nden >= nl + nq.
   OUTPUT
     smult - length nden vector
   ************************************************************ */
void adenscale(double *smult, const double *detd, const int *dencols,
               const int *q, const int *blkend, const int nq, const int nden)
{
  int j,k;
  double detdk;
/* ------------------------------------------------------------
   LORENTZ norm-bound, detd(q(k)) while dencols(j)<blkend(k)
   ------------------------------------------------------------ */
  j= 0;
  for(k = 0; k < nq; k++){
    detdk = detd[q[k]];
    while(j < nden){
      if(dencols[j] >= blkend[k])
        break;
      smult[j++] = detdk;
    }
  }
}

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
  const mxArray *MY_FIELD;
  int i, j, nden, nl, nq, lorN;
  int *q, *dencols, *blkend;
  const double *qPr, *dencolsPr, *detd, *blkstartPr;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  if(nrhs < NPARIN)
    mexErrMsgTxt("adenscale requires more input arguments");
  if(nlhs > NPAROUT)
    mexErrMsgTxt("adenscale produces less output arguments");
/* ------------------------------------------------------------
   DISASSEMBLE dense structure: dense.{l,cols,q}
   ------------------------------------------------------------ */
  if(!mxIsStruct(DENSE_IN))
    mexErrMsgTxt("dense should be a structure.");
  if( (MY_FIELD = mxGetField(DENSE_IN,0,"l")) == NULL)         /* dense.l */
    mexErrMsgTxt("Missing field dense.l.");
  nl = mxGetScalar(MY_FIELD);                           /* double to int */
  if( (MY_FIELD = mxGetField(DENSE_IN,0,"q")) == NULL)         /* dense.q */
    mexErrMsgTxt("Missing field dense.q.");
  nq = mxGetM(MY_FIELD) * mxGetN(MY_FIELD);
  qPr = mxGetPr(MY_FIELD);
  if( (MY_FIELD = mxGetField(DENSE_IN,0,"cols")) == NULL)      /* dense.cols */
    mexErrMsgTxt("Missing field dense.cols.");
  nden = mxGetM(MY_FIELD) * mxGetN(MY_FIELD) - (nq+nl);  /* jump to q-norm */
  if(nden < 0)
    mexErrMsgTxt("dense.cols size mismatch.");
  dencolsPr = mxGetPr(MY_FIELD) + (nq+nl);               /* jump to q-norm */
/* ------------------------------------------------------------
   Disassemble structure d.{det}
   ------------------------------------------------------------ */
  if(!mxIsStruct(D_IN))
    mexErrMsgTxt("d should be a structure.");
  if( (MY_FIELD = mxGetField(D_IN,0,"det")) == NULL)       /* d.det */
    mexErrMsgTxt("Missing field d.det.");
  detd = mxGetPr(MY_FIELD);
  lorN = mxGetM(MY_FIELD) * mxGetN(MY_FIELD);
/* ------------------------------------------------------------
   Get INPUTS blkstart
   ------------------------------------------------------------ */
  blkstartPr = mxGetPr(BLKSTART_IN);
  if(lorN +1 != mxGetM(BLKSTART_IN) * mxGetN(BLKSTART_IN))
    mexErrMsgTxt("blkstart size mismatch");
/* ------------------------------------------------------------
   Create working arrays q(nq), blkend(nq), dencols(nden - (nl+nq))
   ------------------------------------------------------------ */
  q       = (int *) mxCalloc(MAX(nq,1), sizeof(int));
  blkend  = (int *) mxCalloc(MAX(nq,1), sizeof(int));
  dencols = (int *) mxCalloc(MAX(1,nden), sizeof(int));
/* ------------------------------------------------------------
   Convert to integer C-style: q, blkend, dencols.
   ------------------------------------------------------------ */
  for(i = 0; i < nq; i++){
    j = qPr[i];
    q[i] = --j;
  }
  for(i = 0; i < nq; i++){
    j = blkstartPr[q[i] + 1];        /* F-double to C-int */
    blkend[i] = --j;
  }
  for(i = 0; i < nden; i++){
    j = dencolsPr[i];
    dencols[i] = --j;
  }
/* ------------------------------------------------------------
   Create output: smult(nden,1)
   ------------------------------------------------------------ */
  SMULT_OUT = mxCreateDoubleMatrix(nden, 1, mxREAL);
  adenscale(mxGetPr(SMULT_OUT),detd, dencols,q, blkend, nq,nden);
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
  mxFree(dencols);
  mxFree(blkend);
  mxFree(q);
}
