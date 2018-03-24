/*
 Lden = dpr1order(LAD,qloc)
 Creates structure Lden.{LAD,rowperm,colperm,xsuper}

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

*/

#include <string.h>
#include "mex.h"
#include "blksdp.h"

#define LDEN_OUT  plhs[0]
#define NPAROUT 1

#define LAD_IN    prhs[0]
#define QLOC_IN   prhs[1]
#define NPARIN 2

/* ************************************************************
   PROCEDURE dodiscard - discard selected subscripts from selected
      columns.
   INPUT
     xjc - column starts (can be more than m)
     perm - Length m, selection of columns in At.
     discard - discard[xir[inz]] is 1 iff subscript xir[inz] should be
       discarded.
     m - length of perm and xnnz (number of remaining columns).
   UPDATED
     xir - Input: subscript columns, starting at Atjc[perm[j]], of length
       xnnz[j] <= Atjc[perm[j]+1]-Atjc[perm[j]].
       Output: xir(!discard(xir)) adjusted from start of columns
     xnnz - The length of column j, updated after discarding. Length m.
   ************************************************************ */
void dodiscard(int *xir,int *xnnz, const int *xjc,
              const int *perm,const char *discard, const int m)
{
  int i,j,nnzj, inz,jnz;
  int *xj;

  for(j = 0; j < m; j++){
    xj = xir + xjc[perm[j]];
    nnzj = xnnz[j];
    for(inz = 0, jnz = 0; inz < nnzj; inz++)
      if(!discard[(i = xj[inz])])
        xj[jnz++] = i;
    xnnz[j] = jnz;
  }
}

/* ************************************************************
   PROCEDURE dpr1order - Symbolic factorization of dense columns
     (with diag+rank-1 factors) with greedy least-nnz ordering
   INPUT
     m - number of rows (nodes) in x
     n - number of (dense) columns in x
     nqden - number of removed Lorentz blocks. These blocks give an
       additinal column in x.
     xjc - length n+nqden+1, start of each column in x.
     qloc - length nqden, start of removed Lorentz blocks in x(:,1:n).
       We'll let x(:,qloc) = x(:,n+1;n+nqden), and order on the resulting
       first n columns.
   UPDATED
     xir - On input, the row indices of nonzeros in x. Output undefined.
   OUTPUT
     rowperm - Length m. Ordering of the nodes, so that the
       k-th supernode consists of rowperm(xsuper[k]:xsuper[k+1]-1). These
       are the additional nonzeros that factor colperm[k] fans in.
       The entries rowperm(xsuper[n]:m-1) will stay diagonal.
     colperm - Length n. Ordering of the dense column factors, by choosing
       the (first) one with the smallest number of "new" nonzeros.
     xsuper  - Length n+1. Supernodal partition. At stage k, we have
       nonzeros rowperm(0:xsuper[k+1]-1).
   WORK
     xnnz - length n integer working array.
     discard - length m char working array.
   ************************************************************ */
void dpr1order(int *rowperm, int *colperm, int *xsuper, const int *xjc,
               int *xir, const int *qloc, const int m, const int n,
               const int nqden, int *xnnz, char *discard)
{
  int *xk;
  const int *xjck;
  int i,j,k,jnz,pivk, nnzk, permk;
/* ------------------------------------------------------------
   Initialize: colperm = 0:k-1, xnnz(j) = nnz(x(:,j)), discard = all-0,
   xsuper[0]=0.
   ------------------------------------------------------------ */
  for(k = 0; k < n; k++)
    colperm[k] = k;
  for(k = 0; k < n; k++)
    xnnz[k] = xjc[k+1] - xjc[k];
  for(k = 0; k < m; k++)
    discard[k] = 0;
  xsuper[0] = 0;
/* ------------------------------------------------------------
   Permute the last nqden columns with those at qloc
   ------------------------------------------------------------ */
  for(k = 0; k < nqden; k++)
    colperm[qloc[k]] = n+k;
  xjck = xjc + n;                 /* temply point to cols x(:,n+1:n+nqden) */
  for(k = 0; k < nqden; k++){
    j = qloc[k];
    xnnz[j] = xjck[k+1] - xjck[k];
  }
/* ------------------------------------------------------------
   At each main iterate k = 0:n-1, we have:
   xnnz(k:n-1) lists remaining nnz's for remaining cols
   colperm(k:n-1) gives the original column numbers.
   discard(i) = 1 iff row subscript i is in rowperm(0:xsuper[k]-1).
   ------------------------------------------------------------ */
  for(k = 0; k < n; k++){
/* ------------------------------------------------------------
   Pivot on column with smallest number of remaining nonzeros
   ------------------------------------------------------------ */
    nnzk = xnnz[k];
    pivk = k;
    for(j = k+1; j < n; j++)
      if(xnnz[j] < nnzk){
        nnzk = xnnz[j];
        pivk = j;
      }
    permk = colperm[pivk];
    colperm[pivk] = colperm[k];
    colperm[k] = permk;
    xnnz[pivk] = xnnz[k];
/* ------------------------------------------------------------
   Write nonzeros of column "permk" to rowperm and xsuper.
   Then mark them for discard.
   ------------------------------------------------------------ */
    xk = xir + xjc[permk];
    memcpy(rowperm + xsuper[k], xk, nnzk * sizeof(int));
    xsuper[k+1] = xsuper[k] + nnzk;                 /* update xsuper */
    for(i = 0; i < nnzk; i++)
      discard[xk[i]] = 1;
/* ------------------------------------------------------------
   Discard current supernode from remaining columns
   ------------------------------------------------------------ */
    dodiscard(xir,xnnz+(k+1), xjc, colperm+(k+1), discard, n-(k+1));
  }
/* ------------------------------------------------------------
   Write non-discarded row-numbers to bottom of rowperm
   ------------------------------------------------------------ */
  k = 0;
  for(jnz = xsuper[n]; jnz < m; jnz ++){
    for(; discard[k]; k++);
    rowperm[jnz] = k++;                   /* non-discarded node */
  }
  mxAssert(k <= m,"");
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
  const int *xjc;
  int *xir, *rowperm,*colperm,*xsuper, *qloc, *iwork;
  char *cwork;
  double *rowpermPr, *colpermPr, *xsuperPr;
  const double *qlocPr;
  int i,j,m,n,nqden;
  const char *LdenFieldnames[] = {"LAD","rowperm","colperm","xsuper"};
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  if(nrhs < NPARIN)
    mexErrMsgTxt("dpr1order requires more input arguments");
  if(nlhs > NPAROUT)
    mexErrMsgTxt("dpr1order produces less output arguments");
/* ------------------------------------------------------------
   Get column pointers xjc of m x n sparse input matrix X
   ------------------------------------------------------------ */
  if(!mxIsSparse(LAD_IN))
    mexErrMsgTxt("X must be sparse");
  m = mxGetM(LAD_IN);
  n = mxGetN(LAD_IN);
  xjc = mxGetJc(LAD_IN);
/* ------------------------------------------------------------
   Get input qloc: position of Lorentz "d1" columns
   ------------------------------------------------------------ */
  nqden = mxGetM(QLOC_IN) * mxGetN(QLOC_IN);
  qlocPr = mxGetPr(QLOC_IN);
/* ------------------------------------------------------------
   Allocate integer working arrays:
   xir(xjc[n]), rowperm(m), colperm(n-nqden), xsuper(n+1-nqden), qloc(nqden)
   iwork(n-nqden), cwork(m).
   ------------------------------------------------------------ */
  xir     = (int *) mxCalloc(MAX(1,xjc[n]), sizeof(int));
  rowperm = (int *) mxCalloc(MAX(1,m), sizeof(int));
  colperm = (int *) mxCalloc(MAX(1,n-nqden), sizeof(int));
  xsuper  = (int *) mxCalloc(n+1 - nqden, sizeof(int));
  qloc    = (int *) mxCalloc(MAX(nqden,1), sizeof(int));
  iwork   = (int *) mxCalloc(MAX(1,n-nqden), sizeof(int));
  cwork   = (char *) mxCalloc(MAX(1,m), sizeof(char));
/* ------------------------------------------------------------
   Convert qlocPr to int C-style
   ------------------------------------------------------------ */
  for(i = 0; i < nqden; i++){
    j = qlocPr[i];
    qloc[i] = --j;
  }
/* ------------------------------------------------------------
   Initialize xir to subscript structure of X.
   ------------------------------------------------------------ */
  memcpy(xir, mxGetIr(LAD_IN), xjc[n] * sizeof(int));
/* ------------------------------------------------------------
   Then main job is done here:
   ------------------------------------------------------------ */
  dpr1order(rowperm, colperm, xsuper, xjc, xir, qloc, m,n-nqden,nqden,
            iwork,cwork);
/* ------------------------------------------------------------
   Outputs Lden.(LAD, rowperm, colperm, xsuper)
   ------------------------------------------------------------ */
  LDEN_OUT = mxCreateStructMatrix(1, 1, NLDEN_FIELDS, LdenFieldnames);
  LDEN_FIELD = mxDuplicateArray(LAD_IN);               /* LAD */
  mxSetField(LDEN_OUT,0,"LAD",LDEN_FIELD);
  LDEN_FIELD = mxCreateDoubleMatrix(m, 1, mxREAL);     /* rowperm */
  mxSetField(LDEN_OUT,0,"rowperm",LDEN_FIELD);
  rowpermPr = mxGetPr(LDEN_FIELD);
  LDEN_FIELD = mxCreateDoubleMatrix(n, 1, mxREAL);     /* colperm */
  mxSetField(LDEN_OUT,0,"colperm",LDEN_FIELD);
  colpermPr = mxGetPr(LDEN_FIELD);
  LDEN_FIELD  = mxCreateDoubleMatrix(n-nqden+1, 1, mxREAL);  /* xsuper */
  mxSetField(LDEN_OUT,0,"xsuper",LDEN_FIELD);
  xsuperPr = mxGetPr(LDEN_FIELD);
  for(i = 0; i < m; i++)
    rowpermPr[i] = rowperm[i] + 1.0;
  for(i = 0; i <= n-nqden; i++)
    xsuperPr[i] = xsuper[i] + 1.0;
/* ------------------------------------------------------------
   After the rank-1-add columns, we still have the rank-1-subtract
   columns at qloc. Thus, colpermPr = [1+colperm, qlocPr].
   ------------------------------------------------------------ */
  for(i = 0; i < n-nqden; i++)
    colpermPr[i] = colperm[i] + 1.0;
  memcpy(colpermPr + i, qlocPr, nqden * sizeof(double));
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
  mxFree(cwork);
  mxFree(iwork);
  mxFree(qloc);
  mxFree(xsuper);
  mxFree(colperm);
  mxFree(rowperm);
  mxFree(xir);
}
