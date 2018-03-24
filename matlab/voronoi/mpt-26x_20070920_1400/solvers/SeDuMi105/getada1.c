/* ************************************************************
%                                 ADA = getada1(ADA, A,Ajc2,perm, d, blkstart)
% GETADA1  Compute ADA(i,j) = (D(d^2; LP,Lorentz)*A.t(:,i))' *A.t(:,j),
%   and exploit sparsity as much as possible.
%   Ajc2 points just beyond LP/Lorentz nonzeros for each column
%   blkstart = K.qblkstart partitions into Lorentz blocks.
%
%   IMPORTANT 1: only LP and sparse Lorentz part. PSD part ignored altogether.
%   For Lorentz, it uses only det(dk) * ai[k]'*aj[k].
%   IMPORTANT 2: Computes ADA only on triu(ADA(Aord.lqperm,Aord.lqperm)).
%     Remaining entries are set to 0. (CAUTION: sparse(ADA) will therefore
%     destroy the sparsity structure !).
%
% SEE ALSO sedumi, getada2, getada3
% ******************** INTERNAL FUNCTION OF SEDUMI ********************
   function ADA = getada1(ADA, A,Ajc2,perm, d, blkstart)

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

   ************************************************************ */

#include <string.h>
#include <math.h>
#include "mex.h"
#include "blksdp.h"

#define ADA_OUT plhs[0]
#define NPAROUT 1

#define ADA_IN prhs[0]       /* sparsity struct ADA */
#define AT_IN prhs[1]        /* N x m sparse At */
#define AJC2_IN prhs[2]      /* End of LP/Lorentz columns in At */
#define PERM_IN prhs[3]
#define D_IN  prhs[4]         /* scaling vector */
#define BLKSTART_IN  prhs[5]
#define NPARIN 6

/* ************************************************************
   PROCEDURE: getada1
   INPUT
     ada.{jc,ir} - sparsity structure of ada.
     At - sparse N x m matrix.
     d - blkstart[0] (=K.l) vector containing x./z for LP-part.
     ddet - length nblk-1 (=|K.q|) vector containing d.det = (det(dk))_k for
       each Lorentz block k.
     Ajc1 - m int-array, Ajc1 points to start of PSD nz's in At,
       and hence just beyond the LP/Lorentz part.
     blkstart - length nblk+1, cumsum([K.l,|K.q|, K.q-1]), is
         K.blkstart(1:2+length(K.q)).
     perm, invperm - length(m) array, ordering in which ADA should be computed,
       and its inverse. We compute in order triu(ADA(perm,perm)), but store
       at original places. OPTIMAL PERM: sort(Ajc1-At.jc, inc), i.e. start
       with sparsest.
     m  - order of ADA, number of constraints.
     nblk - 1+length(K.q)
   OUTPUT
     ada.pr - ada(i,j) = ai'*D(d^2)*aj. ONLY triu(ADA(perm,perm)) is
        affected. (So caller typically should initialize to all-0.)
   WORKING ARRAYS
     fwork - work vector, size 2*blkstart[nblk].
   ************************************************************ */
void getada1(jcir ada, jcir At, const double *d, const double *ddet,
             const int *Ajc1, const int *blkstart,
             const int *perm, const int *invperm,
             const int m, const int nblk, double *fwork)
{
  int i,j,k, knz,inz, permj;
  double *daj, *dsqr;
  double adaij, detk;
/* ------------------------------------------------------------
   Partition working arrays
   double:  dsqr(lend),  daj(lend), where lend = K.l + sum(K.q.^2).
   ------------------------------------------------------------ */
  daj  = fwork;                                     /* lend */
  dsqr = daj + blkstart[nblk];                      /* lend */
/* ------------------------------------------------------------
   Init daj = all-0 (for LP+Lorentz)
   ------------------------------------------------------------ */
  fzeros(daj, blkstart[nblk]);
/* ------------------------------------------------------------
   Init dsqr = [d.l; -d.det; kron(d.det, all-1)]
   ------------------------------------------------------------ */
  memcpy(dsqr, d, blkstart[0] * sizeof(double));                /* LP */
  memcpy(dsqr + blkstart[0],ddet,(blkstart[1]-blkstart[0]) * sizeof(double));
  for(inz = blkstart[0]; inz < blkstart[1]; inz++)
    dsqr[inz] *= -1;                     /* Lorentz trace */
  ddet -= 2;
  for(k = 2; k <= nblk; k++){            /* Lorentz norm-bound: */
    detk = ddet[k];
    while(inz < blkstart[k])
      dsqr[inz++] = detk;                /* detk * all-1 */
  }
/* ============================================================
   MAIN getada LOOP: loop over nodes perm(0:m-1)
   ============================================================ */
  for(j = 0; j < m; j++){
    permj = perm[j];
    if((inz = At.jc[permj]) < Ajc1[permj]){     /* if any nonzeros */
/* ------------------------------------------------------------
   Compute daj = dsqr .* aj.
   ------------------------------------------------------------ */
      for(; inz < Ajc1[permj]; inz++){
        i = At.ir[inz];
        daj[i] = dsqr[i] * At.pr[inz];
      }
/* ------------------------------------------------------------
   For all i with invpermi < j:
   ada_ij = a_i'*daj.
   ------------------------------------------------------------ */
      for(inz = ada.jc[permj]; inz < ada.jc[permj+1]; inz++){
        i = ada.ir[inz];
        if(invperm[i] <= j){
          for(adaij = 0.0, knz = At.jc[i]; knz < Ajc1[i]; knz++)
            adaij +=  At.pr[knz] * daj[At.ir[knz]];
          ada.pr[inz] = adaij;
        }
      }
/* ------------------------------------------------------------
   Re-initialize daj = 0.
   ------------------------------------------------------------ */
      for(i = At.jc[permj]; i < Ajc1[permj]; i++)   /* LP + Lorentz */
        daj[At.ir[i]] = 0.0;
    } /* ~isempty(At(:,j)) */
  } /* j = 0:m-1 */
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
  int nblk, m, i, j;
  const double *d, *ddet, *permPr, *Ajc2Pr, *blkstartPr;
  double *fwork;
  int *blkstart, *iwork, *Ajc2, *perm, *invperm;
  jcir At, ada;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  if(nrhs < NPARIN)
    mexErrMsgTxt("getADA requires more input arguments.");
  if(nlhs > NPAROUT)
    mexErrMsgTxt("getADA produces less output arguments.");
/* ------------------------------------------------------------
   Get INPUTS blkstart, At, Ajc2, perm.
   ------------------------------------------------------------ */
  blkstartPr = mxGetPr(BLKSTART_IN);
  nblk = mxGetM(BLKSTART_IN) * mxGetN(BLKSTART_IN);   /* is |K.q| + 1 */
  if(nblk < 1)
    mexErrMsgTxt("Size mismatch blkstart.");
  m = mxGetN(AT_IN);
  if(!mxIsSparse(AT_IN))
    mexErrMsgTxt("At should be sparse.");
  At.pr = mxGetPr(AT_IN);
  At.jc = mxGetJc(AT_IN);
  At.ir = mxGetIr(AT_IN);
  Ajc2Pr = mxGetPr(AJC2_IN);
  if(mxGetM(AJC2_IN) * mxGetN(AJC2_IN) != m)
    mexErrMsgTxt("Size mismatch Ajc2.");
  if(mxGetM(PERM_IN) * mxGetN(PERM_IN) != m)
    mexErrMsgTxt("Size mismatch perm.");
  permPr = mxGetPr(PERM_IN);
/* ------------------------------------------------------------
   Allocate working array blkstart(nblk+1).
   ------------------------------------------------------------ */
  blkstart = (int *) mxCalloc(nblk + 1, sizeof(int));
/* ------------------------------------------------------------
   Translate blkstart from Fortran-double to C-int
   ------------------------------------------------------------ */
  for(i = 0; i < nblk; i++){                         /* to integers */
    j = blkstartPr[i];
    blkstart[i+1] = --j;
  }
  if(mxGetM(AT_IN) < blkstart[nblk])
    mexErrMsgTxt("Size mismatch At");
/* ------------------------------------------------------------
   Get SCALING VECTOR: d.{l,det}, and check its size with blkstart.
   ------------------------------------------------------------ */
  if(!mxIsStruct(D_IN))                                       /* d */
    mexErrMsgTxt("Parameter `d' should be a structure.");
  if( (MY_FIELD = mxGetField(D_IN,0,"l")) == NULL)            /* d.l */
    mexErrMsgTxt("Field d.l missing.");
  blkstart[0] = mxGetM(MY_FIELD) * mxGetN(MY_FIELD);
  d = mxGetPr(MY_FIELD);
  if( (MY_FIELD = mxGetField(D_IN,0,"det")) == NULL)            /* d.det */
    mexErrMsgTxt("Field d.det missing.");
  if(mxGetM(MY_FIELD) * mxGetN(MY_FIELD) != blkstart[1] - blkstart[0])
    mexErrMsgTxt("Size d.det mismatch");
  ddet = mxGetPr(MY_FIELD);
/* ------------------------------------------------------------
   Allocate output matrix ADA with sparsity structure of ADA_IN:
   ------------------------------------------------------------ */
  if(mxGetM(ADA_IN) != m || mxGetN(ADA_IN) != m)
    mexErrMsgTxt("Size mismatch ADA.");
  if(!mxIsSparse(ADA_IN))
    mexErrMsgTxt("ADA should be sparse.");
  ada.jc = mxGetJc(ADA_IN);
  ada.ir = mxGetIr(ADA_IN);
  ADA_OUT = mxCreateSparse(m,m, ada.jc[m],mxREAL);  /* ADA = sparse(ADA_IN) */
  ada.pr = mxGetPr(ADA_OUT);                        /* initialized to all-0 */
  memcpy(mxGetJc(ADA_OUT), ada.jc, (m+1) * sizeof(int));
  memcpy(mxGetIr(ADA_OUT), ada.ir, ada.jc[m] * sizeof(int));
/* ------------------------------------------------------------
   ALLOCATE working arrays:
   iwork(3*m) = [Ajc2(m) perm(m), invperm(m)].
   fwork[2 * blkstart[nblk]]
   ------------------------------------------------------------ */
  iwork = (int *) mxCalloc(MAX(3 * m,1), sizeof(int));
  Ajc2 = iwork;
  perm = iwork + m;
  invperm = perm + m;
  fwork  = (double *) mxCalloc(MAX(2 * blkstart[nblk],1), sizeof(double));
/* ------------------------------------------------------------
   perm and Ajc2 to integer C-style
   ------------------------------------------------------------ */
  for(i = 0; i < m; i++){
    j = permPr[i];
    perm[i] = --j;
  }
  for(i = 0; i < m; i++)
    Ajc2[i] = Ajc2Pr[i];
/* ------------------------------------------------------------
   Let invperm(perm) = 0:m-1.
   ------------------------------------------------------------ */
  for(i = 0; i < m; i++)
    invperm[perm[i]] = i;
/* ------------------------------------------------------------
   ACTUAL COMPUTATION: handle constraint aj=At(:,perm(j)), j=0:m-1.
   ------------------------------------------------------------ */
  getada1(ada, At, d, ddet, Ajc2, blkstart, perm, invperm, m, nblk, fwork);
/* ------------------------------------------------------------
   RELEASE WORKING ARRAYS.
   ------------------------------------------------------------ */
  mxFree(fwork);
  mxFree(iwork);
  mxFree(blkstart);
}
