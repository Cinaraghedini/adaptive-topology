/* ************************************************************
function [perm, dz] = dzstruct(At,Ajc1,ifirst)
   
    This file is part of SeDuMi 1.04a
    Copyright (C) 2000 Jos F. Sturm
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
#include "mex.h"
#include "blksdp.h"

#define PERM_OUT myplhs[0]
#define DZ_OUT myplhs[1]
#define NPAROUT 2

#define AT_IN prhs[0]
#define AJC1_IN prhs[1]
#define IFIRST_IN prhs[2]
#define NPARIN 3

/* ************************************************************
   PROCEDURE spPartTransp - Let (Ajc, Air) = At(first:lenfull-1,:)'.
   INPUT
     Atir - length Atjc2[m-1] subscript array.
     Atjc1 - length m array, pointing to 1st index in At(:,j) that is
       in range first:lenfull-1.
     Atjc2 - length m array, pointing to column end.
     first, lenfull - index range we're interested in. (lenud:=lenfull-first)
     m - number of constraints (rows in At).
   OUTPUT
     Ajc - length 1+lenud int-array, column pointers
       of A := At(first:lenfull-1,:)'.
     Air - subscript array of output matrix A. length is:
       sum(Atjc(2:m) - Atjcs(1:m)).
   WORK
     iwork - length lenud (=lenfull-first) integer array.
   ************************************************************ */
void spPartTransp(int *Air, int *Ajc,
                  const int *Atir, const int *Atjc1, const int *Atjc2,
                  const int first, const int lenfull, const int m,
                  int *iwork)
{
  int i,j,inz;
  int *inxnz;
/* ------------------------------------------------------------
   INIT: Make indices Ajc[first:lenfull] valid. Then set Ajc(:)=all-0.
   ------------------------------------------------------------ */
  Ajc -= first;
  for(i = first; i <= lenfull; i++)
    Ajc[i] = 0;
/* ------------------------------------------------------------
   For each column j=0:m-1, in each nz PSD entry Atj(i): ++inxnz[i],
   where inxnz := Ajc+1. (we use Ajc[first+1:lenfull])
   ------------------------------------------------------------ */
  inxnz = Ajc + 1;
  for(j = 0; j < m; j++)
    for(inz = Atjc1[j]; inz < Atjc2[j]; inz++)
      ++inxnz[Atir[inz]];
/* ------------------------------------------------------------
   cumsum(inxnz). Note that Ajc[0] =0 already.
   ------------------------------------------------------------ */
  for(inz = 0, i = first; i < lenfull; i++)
    inxnz[i] += inxnz[i-1];
/* ------------------------------------------------------------
   Write the subscripts of A:= At(first:lenfull-1,:)'.
   ------------------------------------------------------------ */
  memcpy(iwork, Ajc + first, (lenfull - first) * sizeof(int));
  iwork -= first;          /* points to next avl. entry in A(:,i) */
  for(j = 0; j < m; j++)
    for(inz = Atjc1[j]; inz < Atjc2[j]; inz++){
      i = Atir[inz];
      Air[iwork[i]++] = j;         /* as (j,i) entry */
    }
}

/* ************************************************************
   PROCEDURE dzstruct - Greedy ordering of PSD-part of columns in At,
       starting with least number of (PSD)-nonzeros. Good fot getada3.
       Produces N x m sparse matrix dz, having at most lenud nonzeros.
       The last columns correspond to the "deg" columns, and are not
       reordered. The first m-ndeg are reordered, see "perm".
   INPUT
     Atjc1 - length m, start of At(first:lenful,j) for each j=1:m.
     Atjc2 - column end pointers, length m.
     Atir - row-subscripts of At.
     Ajc, Air - sparsity structure of A := At(first:lenful,:)', is m x lenud.
     m    - number of columns in At.
     first - first PSD subscipt.
   OUTPUT
     perm - length m array. perm(1:m-ndeg) is the greedy ordering,
       starting from sparsest. perm(m-ndeg:m) = deg.
     dzjc - length m+1, column pointers into dzir.
     dzir - length dzjc[m] <= lenud. jth column lists additional
       subscripts to dz(:,0:j-1).
   WORK
     iwork - (length m) Remaining length of PSD part
       in each column j=1:m.
     discard - length lenud of Booleans. 1 means index already in dz.
   ************************************************************ */
void dzstruct(int *perm, int *dzjc, int *dzir,
              const int *Atjc1,const int *Atjc2,const int *Atir,
              const int *Ajc,const int *Air,
              const int m, const int first, const int lenud,
              int *iwork, char *discard)
{
  int kmin,lenmin,i,j,k, inz,jnz, permk;
/* ------------------------------------------------------------
   initialize: dzjc[0]=0, nperm = m - ndeg. Let Ajc-=first, so
   that Ajc[first:lenfull] is valid. Similar for discard.
   Initialize discard to all-0 ("False").
   ------------------------------------------------------------ */
  dzjc[0] = 0;
  Ajc -= first;
  memset(discard, '\0', lenud);            /* all-0 */
  discard -= first;
/* ------------------------------------------------------------
   Let perm(1:m) be 1:m
   ------------------------------------------------------------ */
  for(i = 0; i < m; i++)
    perm[i] = i;
/* ------------------------------------------------------------
   Set iwork = Atjc2(1:m)-Atjc1(1:m).
   The number of (not yet discarded) nz-PSD indices per constraint.
   ------------------------------------------------------------ */
  for(k = 0; k < m; k++){
    iwork[k] = Atjc2[k] - Atjc1[k];
  }
/* ------------------------------------------------------------
   In iterate k=0:m-1, pivot on constraint with length iwork[k] minimal.
   ------------------------------------------------------------ */
  for(k = 0; k < m; k++){
    kmin = k;
    lenmin = iwork[perm[k]];
    for(j = k+1; j < m; j++)
      if(iwork[perm[j]] < lenmin){
        lenmin = iwork[perm[j]];
        kmin = j;
      }
    mxAssert(lenmin >= 0,"");
    permk = perm[kmin];                  /* make pivot in perm */
    perm[kmin] = perm[k];
    perm[k] = permk;
/* ------------------------------------------------------------
   Write the (additional) subscripts of the pivot column into
   k-th column of dz, i.e. dz(:,k) = At(:,permk).
   ------------------------------------------------------------ */
    jnz = dzjc[k];
    for(inz = Atjc1[permk]; inz < Atjc2[permk]; inz++){
      i = Atir[inz];
      if(!discard[i]){
        discard[i] = 1;
        dzir[jnz++] = i;
      }
    }
    mxAssert(jnz == dzjc[k] + lenmin,"");
/* ------------------------------------------------------------
   Discard dz(:,k)'s subscripts: adjust the constraint lengths
   where applicable.
   ------------------------------------------------------------ */
    dzjc[k+1] = jnz;
    for(jnz = dzjc[k]; jnz < dzjc[k+1]; jnz++){
      i = dzir[jnz];                         /* i is discarded subscript */
      for(inz = Ajc[i]; inz < Ajc[i+1]; inz++){
        j = Air[inz];
        iwork[j]--;                   /* subscript discard here */
      }
    }
  }
}


/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
      [perm,dz] = dzstruct(At,Ajc1,ifirst)   
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
  mxArray *myplhs[NPAROUT];
  int i,j, m, firstPSD, lenud, iwsiz, lenfull, maxnnz;
  int *iwork, *Atjc1, *Ajc, *Air, *perm;
  const int *Atjc2, *Atir;
  double *permPr;
  const double *Ajc1Pr;
  char *cwork;
  jcir dz;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  if(nrhs < NPARIN)
    mexErrMsgTxt("dzstruct requires more input arguments.");
  if(nlhs > NPAROUT)
    mexErrMsgTxt("dzstruct produces less output arguments.");
/* --------------------------------------------------
   GET STATISTICS:
   -------------------------------------------------- */
  firstPSD = mxGetScalar(IFIRST_IN);        /* double to int */
  --firstPSD;                               /* Fortron to C */
  lenfull = mxGetM(AT_IN);
  lenud = lenfull - firstPSD;
/* --------------------------------------------------
   Check size At, and get At.
   -------------------------------------------------- */
  if(!mxIsSparse(AT_IN))
    mexErrMsgTxt("At must be a sparse matrix.");
  m = mxGetN(AT_IN);
  Atjc2 = mxGetJc(AT_IN)+1;        /* points to end of constraint */
  Atir = mxGetIr(AT_IN);           /* subscripts */
/* ------------------------------------------------------------
   Get input AJc1
   ------------------------------------------------------------ */
  Ajc1Pr = mxGetPr(AJC1_IN);
  if(mxGetM(AJC1_IN) * mxGetN(AJC1_IN) < m)
    mexErrMsgTxt("Ajc1 size mismatch");
/* ------------------------------------------------------------
   ALLOCATE WORKING arrays:
   int Atjc1(m),  Ajc(1+lenud), Air(maxnnz), perm(m),
     iwork(iwsiz). iwsiz = MAX(lenud,1+m)
   char cwork(lenud).
   ------------------------------------------------------------ */
  Atjc1 = (int *) mxCalloc(MAX(m,1), sizeof(int));
  Ajc = (int *) mxCalloc(1+lenud, sizeof(int));
  perm = (int *) mxCalloc(MAX(1,m), sizeof(int));
  iwsiz = MAX(lenud, 1 + m);
  iwork = (int *) mxCalloc(iwsiz, sizeof(int)); /* iwork(iwsiz) */
  cwork = (char *) mxCalloc(MAX(1,lenud), sizeof(char));
/* ------------------------------------------------------------
   Double to int: Atjc1
   ------------------------------------------------------------ */
  for(i = 0; i < m; i++)
    Atjc1[i] = Ajc1Pr[i];           /* double to int */
/* ------------------------------------------------------------
   Let maxnnz = number of PSD nonzeros = sum(Atjc2-Atjc1)
   ------------------------------------------------------------ */
  maxnnz = 0;
  for(i = 0; i < m; i++)
    maxnnz += Atjc2[i] - Atjc1[i];
/* ------------------------------------------------------------
   ALLOCATE WORKING array int Air(maxnnz)
   ------------------------------------------------------------ */
  Air = (int *) mxCalloc(MAX(1,maxnnz), sizeof(int));
/* ------------------------------------------------------------
   CREATE OUTPUT ARRAYS PERM(m) and DZ=sparse(lenfull,m,lenud).
   ------------------------------------------------------------ */
  PERM_OUT = mxCreateDoubleMatrix(m,1,mxREAL);
  permPr = mxGetPr(PERM_OUT);
  DZ_OUT = mxCreateSparse(lenfull,m,lenud,mxREAL);
  dz.jc = mxGetJc(DZ_OUT);
  dz.ir = mxGetIr(DZ_OUT);
/* ------------------------------------------------------------
   Let (Ajc,Air) := At(first:end,:)', the transpose of PSD-part.
   Uses iwork(lenud)
   ------------------------------------------------------------ */
  spPartTransp(Air,Ajc, Atir,Atjc1,Atjc2, firstPSD,lenfull, m,iwork);
/* ------------------------------------------------------------
   The main job: greedy order of columns of At
   ------------------------------------------------------------ */
  dzstruct(perm, dz.jc,dz.ir, Atjc1,Atjc2,Atir, Ajc,Air,
           m, firstPSD,lenud, iwork, cwork);    /* uses iwork(m+1) */
/* ------------------------------------------------------------
   REALLOC (shrink) dz to dz.jc[m] nonzeros.
   ------------------------------------------------------------ */
  mxAssert(dz.jc[m] <= lenud,"");
  maxnnz = MAX(1,dz.jc[m]);                     /* avoid realloc to 0 */
  if((dz.ir = (int *) mxRealloc(dz.ir, maxnnz * sizeof(int))) == NULL)
    mexErrMsgTxt("Memory allocation error");
  if((dz.pr = (double *) mxRealloc(mxGetPr(DZ_OUT), maxnnz*sizeof(double)))
     == NULL)
    mexErrMsgTxt("Memory allocation error");
  mxSetPr(DZ_OUT,dz.pr);
  mxSetIr(DZ_OUT,dz.ir);
  mxSetNzmax(DZ_OUT,maxnnz);
  for(i = 0; i < maxnnz; i++)
    dz.pr[i] = 1.0;
/* ------------------------------------------------------------
   Convert C-int to Fortran-double
   ------------------------------------------------------------ */
  for(i = 0; i < m; i++)
    permPr[i] = perm[i] + 1.0;
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
  mxFree(cwork);
  mxFree(iwork);
  mxFree(perm);
  mxFree(Air);
  mxFree(Ajc);
  mxFree(Atjc1);
/* ------------------------------------------------------------
   Copy requested output parameters (at least 1), release others.
   ------------------------------------------------------------ */
  i = MAX(nlhs, 1);
  memcpy(plhs,myplhs, i * sizeof(mxArray *));
  for(; i < NPAROUT; i++)
    mxDestroyArray(myplhs[i]);
}
