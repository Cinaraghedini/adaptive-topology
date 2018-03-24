/*
%                            Ad = Adendotd(dense, d, sparAd, Ablk, blkstart)
% ADENDOTD  Computes d[k]'*Aj[k] for Lorentz blocks that are to be factored
%  by dpr1fact.
%
% SEE ALSO sedumi
% **********  INTERNAL FUNCTION OF SEDUMI **********
function Ad = Adendotd(dense, d, sparAd, Ablk, blkstart)

    This file is part of SeDuMi 1.05
    Copyright (C) 2001 Jos F. Sturm
      Dept. Econometrics & O.R., Tilburg University, the Netherlands.
      Supported by the Netherlands Organization for Scientific Research (NWO).
    Affiliation SeDuMi 1.03 and 1.04Beta (2000):
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
#include <string.h>
#include "mex.h"
#include "blksdp.h"

#define AD_OUT plhs[0]
#define NPAROUT 1

#define DENSE_IN prhs[0]
#define D_IN prhs[1]
#define ADOTD_IN prhs[2]
#define ABLK_IN prhs[3]
#define BLKSTART_IN prhs[4]
#define NPARIN 5

/* ************************************************************
   PROCEDURE adendotd
   INPUT
     aden - dense.A(:,dense.l+1:end) sparse m x (nq+nden) matrix,
       with aden.jc[0] possibly nonzero.
     adotd - sparse m x nq array, has ai[k]'*d[k] for k in q, where
       the ai's are the sparse part of the A-matrix. We still need to
       add contribution from dense part, resulting in Ad.
     d1 - length |K.q| vector. We will use d1(q) entries.
     d2 - length firstQ+(sym(K.q)-|K.q|) vector. We use entries d2(dencols),
       where dencols >=firstQ.
     q - length nq array: dense lorentz blocks
     dencols - length nden array: dense lorentz norm-bound columns. These
       are global subscripts, at or beyond firstQ.
     blkend - length nq array, listing 1-beyond-last subscript of Lorentz
       norm bound blocks listed in q.
     nq - number of dense lorentz blocks
     nden - number of dense lorentz norm-bound columns
     fwork - length m vector.
   UPDATED
     ad - sparse m x nq. ad.ir and ad.jc are INPUTS, ad.pr is OUTPUT.
       On output, has (ai[k]+Adeni[k])'*d[k] for k in q.
   ************************************************************ */
void adendotd(jcir ad,jcir adotd,jcir aden,const double *d1,const double *d2,
              const int *q,const int *dencols,
              const int *blkend,const int nq,const int nden, double *fwork)
{
  int inz, i,j,k;
  const int *aden2jc;
  double dj;
/* ------------------------------------------------------------
   Initialize (Lorentz norm-bound part):
   1) aden2jc(0:nden) points to dense columns 
   2) j is next dense column to handle, inz point to next nonzero
   ------------------------------------------------------------ */
  j = 0;
  aden2jc = aden.jc + nq;   /* jump over Lorentz trace columns*/
  inz = aden2jc[j];
  for(k = 0; k < nq; k++){
/* ------------------------------------------------------------
   Set fwork = all-0;
   ------------------------------------------------------------ */
    for(i = ad.jc[k]; i < ad.jc[k+1]; i++)        /* fwork = all-0 */
      fwork[ad.ir[i]] = 0.0;
/* ------------------------------------------------------------
   Let fwork = adotd(:,k)    (Contribution from SPARSE part of A)
   ------------------------------------------------------------ */
    for(i = adotd.jc[k]; i < adotd.jc[k+1]; i++)
      fwork[adotd.ir[i]] = adotd.pr[i];
/* ------------------------------------------------------------
   Let fwork += d1(q(k)) * aden(:,k)   (Contribution Lorentz-trace)
   ------------------------------------------------------------ */
    dj = d1[q[k]];
    for(i = aden.jc[k]; i < aden.jc[k+1]; i++)
      fwork[aden.ir[i]] += dj * aden.pr[i];
/* ------------------------------------------------------------
   Add contribution of dense Lorentz-norm-bound columns, i.e.
   let fwork += sum_j{d2(dencols[j]) * Aden(:,j) | dencols[j]<blkend[k]}
   ------------------------------------------------------------ */
    for(; j < nden; j++){
      if((i = dencols[j]) >= blkend[k])
        break;                            /* Break if beyond block k */
      dj = d2[i];
      for(; inz < aden2jc[j+1]; inz++)
        fwork[aden.ir[inz]] += dj * aden.pr[inz];
    }
/* ------------------------------------------------------------
   Store ad(:,k) = fwork
   ------------------------------------------------------------ */
    for(i = ad.jc[k]; i < ad.jc[k+1]; i++)
      ad.pr[i] = fwork[ad.ir[i]];
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
  int i,j,firstQ, m,nden, nl, nq, lorN;
  int *q, *dencols, *blkend;
  const double *d1, *d2, *qPr, *dencolsPr, *blkstartPr;
  double *fwork;
  jcir ad, aden,adotd;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  if(nrhs < NPARIN)
    mexErrMsgTxt("adendotd requires more input arguments");
  if(nlhs > NPAROUT)
    mexErrMsgTxt("adendotd produces less output arguments");
/* ------------------------------------------------------------
   DISASSEMBLE dense structure: dense.{cols,l,q,A}
   ------------------------------------------------------------ */
  if(!mxIsStruct(DENSE_IN))
    mexErrMsgTxt("dense should be a structure.");
  if( (MY_FIELD = mxGetField(DENSE_IN,0,"l")) == NULL)      /* dense.l */
    mexErrMsgTxt("Missing field dense.l.");
  nl = mxGetScalar(MY_FIELD);                             /* double to int */
  if( (MY_FIELD = mxGetField(DENSE_IN,0,"q")) == NULL)      /* dense.q */
    mexErrMsgTxt("Missing field dense.q.");
  nq = mxGetM(MY_FIELD) * mxGetN(MY_FIELD);
  qPr = mxGetPr(MY_FIELD);
  if( (MY_FIELD = mxGetField(DENSE_IN,0,"cols")) == NULL)      /* dense.cols */
    mexErrMsgTxt("Missing field dense.cols.");
  nden = mxGetM(MY_FIELD) * mxGetN(MY_FIELD) - nl - nq;
  if(nden < 0)
    mexErrMsgTxt("dense.q size mismatch.");
  dencolsPr = mxGetPr(MY_FIELD) + nl + nq;               /* Skip LP and Q-tr*/
  if( (MY_FIELD = mxGetField(DENSE_IN,0,"A")) == NULL)      /* dense.A */
    mexErrMsgTxt("Missing field dense.A.");
  if(!mxIsSparse(MY_FIELD))
    mexErrMsgTxt("dense.A must be sparse");
  m = mxGetM(MY_FIELD); 
  if(mxGetN(MY_FIELD) - nl != nq + nden)
    mexErrMsgTxt("dense.A size mismatch");
  aden.jc = mxGetJc(MY_FIELD) + nl;                         /* Skip LP part */
  aden.ir = mxGetIr(MY_FIELD);
  aden.pr = mxGetPr(MY_FIELD);
/* ------------------------------------------------------------
   DISASSEMBLE d structure: d.{q1,q2}
   ------------------------------------------------------------ */
  if(!mxIsStruct(D_IN))
    mexErrMsgTxt("d should be a structure.");
  if( (MY_FIELD = mxGetField(D_IN,0,"q1")) == NULL)      /* d.q1 */
    mexErrMsgTxt("Missing field d.q1.");
  lorN = mxGetM(MY_FIELD) * mxGetN(MY_FIELD);
  d1 = mxGetPr(MY_FIELD);
  if( (MY_FIELD = mxGetField(D_IN,0,"q2")) == NULL)      /* d.q2 */
    mexErrMsgTxt("Missing field d.q2.");
  d2 = mxGetPr(MY_FIELD);
/* ------------------------------------------------------------
   Get inputs adotd (contains Ad from sparse A in dense.qs blocks),
   blkstart (partitions d2 into Lorentz norm-bound blocks)
   ------------------------------------------------------------ */
  if(!mxIsSparse(ADOTD_IN))                            /* adotd */
    mexErrMsgTxt("sparAD must be sparse");
  if((m != mxGetM(ADOTD_IN) && nq > 0) || nq != mxGetN(ADOTD_IN))
    mexErrMsgTxt("Size mismatch sparAD");
  adotd.jc = mxGetJc(ADOTD_IN);
  adotd.ir = mxGetIr(ADOTD_IN);
  adotd.pr = mxGetPr(ADOTD_IN);
  blkstartPr = mxGetPr(BLKSTART_IN);                  /* blkstart */
  if(lorN +1 != mxGetM(BLKSTART_IN) * mxGetN(BLKSTART_IN))
    mexErrMsgTxt("blkstart size mismatch");
/* ------------------------------------------------------------
   Create working arrays q(nq), dencols(nden), fwork(m),
   blkend(nq)
   ------------------------------------------------------------ */
  q    = (int *) mxCalloc(MAX(1,nq), sizeof(int));
  dencols = (int *) mxCalloc(MAX(1,nden), sizeof(int));
  blkend = (int *) mxCalloc(MAX(1,nq), sizeof(int));
  fwork = (double *) mxCalloc(MAX(m,1), sizeof(double));
/* ------------------------------------------------------------
   Convert to integer C-style; dencols, q, blkstart(q+1)
   ------------------------------------------------------------ */
  for(i = 0; i < nden; i++){
    j = dencolsPr[i];
    dencols[i] = --j;
  }
  for(i = 0; i < nq; i++){
    j = qPr[i];
    q[i] = --j;
  }
/* ------------------------------------------------------------
   Let firstQ point to subscript of 1st Lorentz norm-bound variable
   ------------------------------------------------------------ */
  firstQ = blkstartPr[0];            /* double to int */
  --firstQ;                          /* Fortran to C */
  for(i = 0; i < nq; i++){
    j = blkstartPr[q[i] + 1];        /* F-double to C-int */
    blkend[i] = --j;
  }
/* ------------------------------------------------------------
   Create output: Ad = Ablk
   ------------------------------------------------------------ */
  AD_OUT = mxDuplicateArray(ABLK_IN);              /* Ad = Ablk */
  ad.jc = mxGetJc(AD_OUT);
  ad.ir = mxGetIr(AD_OUT);
  ad.pr = mxGetPr(AD_OUT);
/* ------------------------------------------------------------
   The real job is done here:
   ------------------------------------------------------------ */
  adendotd(ad,adotd,aden,d1,d2 - firstQ,q,dencols,blkend,nq,nden, fwork);
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
  mxFree(fwork);
  mxFree(dencols);
  mxFree(q);
}
