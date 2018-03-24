/*
%                                   Lden = finsymbden(LAD,perm,dz,firstq)
% FINSYMBDEN  Updates perm and dz by inserting the
%  last Lorentz trace columns (last columns of LAD). It creates the fields
%  Lden.sign  - +1 for "normal" columns, -1 for Lorentz trace columns
%  Lden.first - First pivot column that will affect this one
%  NOTE: sign and first correspond to columns in LAD (without perm-reordering).
%
% SEE ALSO incorder
% ******************** INTERNAL FUNCTION OF SEDUMI ********************
function Lden = finsymbden(LAD,perm,dz,firstq)

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

#define LDEN_OUT  plhs[0]
#define NPAROUT 1

#define LAD_IN    prhs[0]
#define PERM_IN   prhs[1]
#define DZ_IN     prhs[2]
#define FIRSTQ_IN prhs[3]
#define NPARIN 4

/* ************************************************************
   PROCEDURE getfirstpiv - Find first affecting pivot on column
     j = 1:n. These column numbers are in NON-PIVOTED ORDER, i.e.
     the order in which they appear in Xjc.
   INPUT
     invperm - length n array, yields position in list of nonzeros "dzir".
     xsuper - length n (though it may have n+1) array, partitioning
       of permuted subscripts, is "dzjc".
     Xjc - length n+1 array
     Xir - length Xjc[n] array
     n - number of columns in X.
   OUTPUT
     firstpiv - length n array
   ************************************************************ */
void getfirstpiv(int *firstpiv, const int *invperm, const int *xsuper,
                 const int *Xjc, const int *Xir, const int n)
{
  int i,j,inz,firstj;
  inz = Xjc[0];                 /* typically inz = 0*/
  for(j = 0; j < n; j++){
/* ------------------------------------------------------------
   Let firstj = min(invperm(find(X(:,j))))
   ------------------------------------------------------------ */
    if(inz < Xjc[j+1]){
      firstj = invperm[Xir[inz]];
      for(++inz; inz < Xjc[j+1]; inz++)
        if((i = invperm[Xir[inz]]) < firstj)
          firstj = i;
/* ------------------------------------------------------------
   First node covering firstj, i.e. xsuper[y] < firstj+1 <= xsuper[y+1],
   with y denoting firstpiv[j].
   ------------------------------------------------------------ */
      firstpiv[j] = 0;             /* search from start */
      intbsearch(firstpiv+j,xsuper+1,n-1,firstj+1);
    }
    else
      firstpiv[j] = n;             /* if all-0 then no affecting pivot */
  }
}

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
#define NLDEN_FIELDS 4
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
  mxArray *LDEN_FIELD;
  int i,inz, j,m,n,nperm, firstQ, lastQ, nnzdz;
  const int *LADjc, *LADir, *dzJc, *dzIr;
  int *invdz,*firstpiv,*perm, *dznewJc;
  double *permPr, *firstPr;
  const char *LdenFieldnames[] = {"LAD","perm","dz","first"};
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  if(nrhs < NPARIN)
    mexErrMsgTxt("finsymbden requires more input arguments");
  if(nlhs > NPAROUT)
    mexErrMsgTxt("finsymbden produces less output arguments");
/* ------------------------------------------------------------
   Get inputs LAD,perm,dz,firstq
   n = total number of dense columns/blocks, i.e. #cols in LAD
   m = number of constraints
   nperm = n - number of removed lorentz trace columns.
   ------------------------------------------------------------ */
  if(!mxIsSparse(LAD_IN))                      /* LAD */
    mexErrMsgTxt("LAD must be sparse");
  m = mxGetM(LAD_IN);
  n = mxGetN(LAD_IN);
  LADjc = mxGetJc(LAD_IN);
  LADir = mxGetIr(LAD_IN);
  permPr = mxGetPr(PERM_IN);                    /* perm */
  nperm = mxGetM(PERM_IN) * mxGetN(PERM_IN);
  dzJc = mxGetJc(DZ_IN);                        /* dz */
  dzIr = mxGetIr(DZ_IN);
  if(mxGetM(DZ_IN) != m || mxGetN(DZ_IN) != nperm)
    mexErrMsgTxt("dz size mismatch");
/* ------------------------------------------------------------
   INPUT firstQ == dense.l+1, points to 1st entry in dense.cols
    dealing with Lorentz-trace entries. Let lastQ point just beyond
    Lorentz trace/block entries, i.e. add n-nperm.
   ------------------------------------------------------------ */
  firstQ = mxGetScalar(FIRSTQ_IN);         /*firstq, F-double to C-int.*/
  --firstQ;
  lastQ = firstQ + n - nperm;
/* ------------------------------------------------------------
   Allocate integer working arrays:
   invdz(m), firstpiv(n), perm(n)
   ------------------------------------------------------------ */
  invdz = (int *) mxCalloc(MAX(1,m), sizeof(int));
  firstpiv = (int *) mxCalloc(MAX(1,n), sizeof(int));
  perm = (int *) mxCalloc(MAX(1,n), sizeof(int));
/* ------------------------------------------------------------
   Allocate OUTPUT int array dznewJc(n+1)
   ------------------------------------------------------------ */
  dznewJc = (int *) mxCalloc(n+1, sizeof(int));
/* ------------------------------------------------------------
   Let invdz(dzIr) = 1:nnz(dz). Note that nnz(dz)<m is the number
   subscripts that are actually in use.
   ------------------------------------------------------------ */
  nnzdz = dzJc[nperm];
  for(i = dzJc[0]; i < nnzdz; i++)
    invdz[dzIr[i]] = i;                 /* dz is m x nperm */
/* ------------------------------------------------------------
   Create new perm and dz-column pointers, to include lorentz trace cols.
   These cols are attached to Lorentz-blocks cols, whose subscripts
   range in firstQ:lastQ-1.
   ------------------------------------------------------------ */
  inz = 0;
  for(i = 0; i < nperm; i++){
    j = permPr[i];
    perm[inz] = --j;
    dznewJc[inz++] = dzJc[i];
    if(j >= firstQ && j < lastQ){
/* ------------------------------------------------------------
   Attach Lorentz trace col. These cols are at nperm:n-1.
   ------------------------------------------------------------ */
      perm[inz] = nperm + j - firstQ;   /* insert associated trace column */
      mxAssert(perm[inz] < n,"");
      dznewJc[inz++] = dzJc[i+1];      /* no extra subscripts->start at end */
    }
  }
  mxAssert(inz == n,"");
  dznewJc[n] = dzJc[nperm];
/* ------------------------------------------------------------
   Compute firstpiv
   ------------------------------------------------------------ */
  getfirstpiv(firstpiv, invdz, dznewJc, LADjc,LADir, n);
/* ------------------------------------------------------------
   Outputs Lden.(LAD, perm, dz, first)
   ------------------------------------------------------------ */
  LDEN_OUT = mxCreateStructMatrix(1, 1, NLDEN_FIELDS, LdenFieldnames);
  LDEN_FIELD = mxDuplicateArray(LAD_IN);               /* LAD */
  mxSetField(LDEN_OUT,0,"LAD",LDEN_FIELD);
  LDEN_FIELD = mxCreateDoubleMatrix(n, 1, mxREAL);     /* perm */
  mxSetField(LDEN_OUT,0,"perm",LDEN_FIELD);
  permPr = mxGetPr(LDEN_FIELD);
  for(i = 0; i < n; i++)
    permPr[i] = perm[i] + 1.0;                       /* C-int to F-double */
  LDEN_FIELD = mxDuplicateArray(DZ_IN);                /* dz */
/* NOTE: here we replace jc by dznewJc */
  mxFree(mxGetJc(LDEN_FIELD));
  mxSetJc(LDEN_FIELD, dznewJc);
  mxSetN(LDEN_FIELD, n);
  mxSetField(LDEN_OUT,0,"dz",LDEN_FIELD);
  LDEN_FIELD  = mxCreateDoubleMatrix(n, 1, mxREAL);  /* first */
  mxSetField(LDEN_OUT,0,"first",LDEN_FIELD);
  firstPr = mxGetPr(LDEN_FIELD);
  for(i = 0; i < n; i++)
    firstPr[i] = firstpiv[i] + 1.0;               /* C-int to F-double */
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
  mxFree(perm);
  mxFree(firstpiv);
  mxFree(invdz);
}
