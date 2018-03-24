/*
 x = symbfwblk(L, b)

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

#define X_OUT  plhs[0]
#define NPAROUT 1

#define L_IN      prhs[0]
#define B_IN      prhs[1]
#define NPARIN 2

/* ************************************************************
   PROCEDURE snodeCompress  - Compressed subscripts based on
                    supernodal partition.
   INPUT
     ljc, lir - uncompressed nz-structure of L
     xsuper - supernodal partition (length nsuper+1).
     nsuper - number of supernodes
   OUTPUT
     xlindx - length nsuper+1, the columns in lindx.
     lindx  - compressed subscript array: L(:,xsuper). Should be allocated
       to ljc[m], so that there are certainly enough entries.
     snode  - length m = xsuper[nsuper+1]. Maps subnode to its supernode.
   ************************************************************ */
void snodeCompress(int *xlindx,int *lindx,int *snode,
                   const int *ljc,const int *lir,const int *xsuper,
                   const int nsuper)
{
  int j, jsup, ix, collen, jcol;
/* ------------------------------------------------------------
   SNODE: map each column to the supernode containing it
   ------------------------------------------------------------ */
  j = xsuper[0];
  for(jsup = 0; jsup < nsuper; jsup++){
    while(j < xsuper[jsup + 1])
      snode[j++] = jsup;
  }
/* ------------------------------------------------------------
   COMPRESS SUBSCRIPTS:
    Let (xlindx,lindx) = ljc(xsuper(:)), i.e store only once
    for each snode, instead of once per column.
   ------------------------------------------------------------ */
  for(ix = 0, jsup = 0; jsup < nsuper; jsup++){
    xlindx[jsup] = ix;
    jcol = xsuper[jsup];
    collen = ljc[jcol+1] - ljc[jcol];
    memcpy(lindx + ix, lir + ljc[jcol], collen * sizeof(int));
    ix += collen;
  }
  xlindx[nsuper] = ix;
}

/* ************************************************************
   PROCEDURE getnzfwlj - find nonzero supernodes in L\e_{xsuper[jsup]}.
   INPUT
     jsup  - starting supernode to process from.
     snode, xsuper - supernodal partition.
         xsuper(nsuper+1): xsuper(j):xsuper(j+1)-1 is jth supernode
         snode(m): j=snode(i) means xsuper(j) <= i < xsuper(j+1).
     xlindx,lindx  - compressed subscript array.
         xlindx(nsuper+1): lindx(xlindx(j):xlindx(j+1)-1) are the subscripts
         for supernode j.
   UPDATED
     processed - Sets processed[jsup] = 1 if jsup is in L\e_{xsuper[jsup]}.
     snodefrom - Lists first relevant subnode of each supernode 
       where processed[jsup]=1.
   REMARK - caller has to set processed[jsup]=1; getnzfwlj only does
     this for the child supernodes.
   ************************************************************ */
void getnzfwlj(int *snodefrom, char *processed, int jsup,
              const int *snode, const int *xsuper,
              const int *xlindx, const int *lindx)
{
  int i,j;

  j = xsuper[jsup+1] - xsuper[jsup];
  while(xlindx[jsup] + j < xlindx[jsup + 1]){
    i = lindx[xlindx[jsup] + j];
    jsup = snode[i];         /* next affected snode */
/* ------------------------------------------------------------
   If jsup has already been processed, then we can stop here, after
   making sure that i >= snodefrom[jsup].
   ------------------------------------------------------------ */
    if(processed[jsup]){
      if(i < snodefrom[jsup])
        snodefrom[jsup] = i;
      break;             /* STOP */
    }
/* ------------------------------------------------------------
   Otherwise, we process and link through to next affected supernode
   ------------------------------------------------------------ */
    processed[jsup] = 1;
    snodefrom[jsup] = i;
    j = xsuper[jsup+1] - xsuper[jsup];
  }
}

/* ************************************************************
   PROCEDURE getnzsuper - Compute sparsity structure of L\b(perm), by
       determining nonzero-supernodes, and starting subnodes within
       them (each supernode is a dense diag block in L).
   INPUT
     bir, bnnz     - bir(bnnz) lists the row-indices of vector b.
     invperm       - int(m) Is s.t. perm[invperm[i]] = i.
     snode, xsuper - supernodal partition.
         xsuper(nsuper+1): xsuper(j):xsuper(j+1)-1 is jth supernode
         snode(m): j=snode(i) means xsuper(j) <= i < xsuper(j+1).
     xlindx,lindx  - compressed subscript array.
         xlindx(nsuper+1): lindx(xlindx(j):xlindx(j+1)-1) are the subscripts
         for supernode j.
   UPDATED
     processed - char(nsuper) array. On input all-0, on output
       processed[jsup] = 1 iff jsup is in nz structure of L\b.
   OUTPUT
     snodefrom - Length nsuper array. Lists first relevant subnode of
       each supernode where processed[jsup]=1.
   ************************************************************ */
void getnzsuper(int *snodefrom, char *processed,
                const int *bir, const int bnnz,
                const int *invperm, const int *snode, const int *xsuper,
                const int *xlindx, const int *lindx)
{
  int inz,i,jsup;
/* ------------------------------------------------------------
   We'll process each supernode jsup = snode[ bir[ inz ] ], to
   find all supernodes in x, L*x = b, and the first relevant
   subnode of each supernode, snodefrom[jsup].
   ------------------------------------------------------------ */
  if(bnnz <= 0)
    return;
  for(inz = 0; inz < bnnz; inz++){
    i = invperm[bir[inz]];            /* We're interested in b(perm) */
    jsup = snode[i];
/* ------------------------------------------------------------
   If jsup has not yet been processed, then find supernodes involved
   in solving L*x = e_i.
   ------------------------------------------------------------ */
    if(!processed[jsup]){
      snodefrom[jsup] = i;
      getnzfwlj(snodefrom,processed,jsup, snode,xsuper,xlindx,lindx);
      processed[jsup] = 1;
    }
/* ------------------------------------------------------------
   Otherwise, we only need to make sure that i >= snodefrom[jsup].
   ------------------------------------------------------------ */
    else if(i < snodefrom[jsup])
      snodefrom[jsup] = i;
  }
}    

/* ************************************************************
   PROCEDURE symbfwmat - Computes nz-structur of x = L\b(perm,:).
   INPUT
     bjc, bir - nz-structure of m x n RHS-matrix b.
     invperm       - int(m) Is s.t. perm[invperm[i]] = i.
     snode, xsuper - supernodal partition.
         xsuper(nsuper+1): xsuper(j):xsuper(j+1)-1 is jth supernode
         snode(m): j=snode(i) means xsuper(j) <= i < xsuper(j+1).
     xlindx,lindx  - compressed subscript array.
         xlindx(nsuper+1): lindx(xlindx(j):xlindx(j+1)-1) are the subscripts
         for supernode j.
     nsuper   - number of supernodes
     m,n      - size(b), m rows, n columns.
   OUTPUT
     xjc     - n+1-array: Start of each column in *pxir.
     *pxir   - length *pmaxnnz array of row-indices.
   UPDATED
     *pmaxnnz - The allocated number of entries in *pxir. Will be changed
        by this function to the exact number needed (but s.t. maxnnz >= 1).
   WORK
     snodefrom - int(nsuper)
     processed - char(nsuper)
   ************************************************************ */
void symbfwmat(int *xjc, int **pxir,int *pmaxnnz,
               const int *bjc, const int *bir,
               const int *invperm, const int *snode, const int *xsuper,
               const int *xlindx, const int *lindx,
               const int nsuper, const int m, const int n,
               int *snodefrom, char *processed)
{
  int i,j,k,inz, maxnnz;
  int *xir;
/* ------------------------------------------------------------
   INIT: processed = 0, xir = *pxir, maxnnz = *pmaxnnz.
   ------------------------------------------------------------ */
  memset(processed, 0, nsuper * sizeof(char));
  xir = *pxir;
  maxnnz = *pmaxnnz;
/* ------------------------------------------------------------
   For each column j, compute nz-structure of L\b(:,j) into xir.
   First make sure that xir has enough (at least m) available entries.
   ------------------------------------------------------------ */
  inz = 0;
  for(j = 0; j < n; j++){
    xjc[j] = inz;
    if(inz + m > maxnnz){
      maxnnz += inz + m;        /* required + old amount */
      xir = (int *) mxRealloc(xir, maxnnz*sizeof(int));
    }
/* ------------------------------------------------------------
   Find all nz-supernodes in L\b(:,j).
   ------------------------------------------------------------ */
    getnzsuper(snodefrom, processed, bir + bjc[j], bjc[j+1]-bjc[j],
               invperm, snode, xsuper, xlindx, lindx);
/* ------------------------------------------------------------
   For each nz-supernode, write the row-indices from "snodefrom" into xir.
   ------------------------------------------------------------ */
    for(k = 0; k < nsuper; k++)
      if(processed[k]){
        processed[k] = 0;
        for(i = snodefrom[k]; i < xsuper[k+1]; i++)
          xir[inz++] = i;
      }
  }
/* ------------------------------------------------------------
   FINALLY: close last column in xir, Realloc xir to the actual
   maxnnz:= xjc[n], and return.
   ------------------------------------------------------------ */
  xjc[n] = inz;
  if(inz < maxnnz){
    maxnnz = MAX(inz,1);           /* avoid realloc to NULL */
    xir = (int *) mxRealloc(xir, maxnnz * sizeof(int));
  }
  *pxir = xir;
  *pmaxnnz = maxnnz;
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
  const mxArray *L_FIELD;
  int maxnnz, i,j, nsuper,m,n;
  const int *ljc,*lir,*bjc,*bir;
  int *xjc,*xir, *snode,*xlindx,*lindx, *iwork,*xsuper, *invperm;
  char *cwork;
  double *xpr;
  const double *permPr, *xsuperPr;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  if(nrhs < NPARIN)
    mexErrMsgTxt("symbfwblk requires more input arguments");
  if(nlhs > NPAROUT)
    mexErrMsgTxt("symbfwblk produces 1 output argument");
/* ------------------------------------------------------------
   Get rhs-input B
   ------------------------------------------------------------ */
  if(!mxIsSparse(B_IN))
    mexErrMsgTxt("B must be sparse");
  m = mxGetM(B_IN);
  n = mxGetN(B_IN);
  bjc = mxGetJc(B_IN);
  bir = mxGetIr(B_IN);
/* ------------------------------------------------------------
   Disassemble block Cholesky structure L
   ------------------------------------------------------------ */
  if(!mxIsStruct(L_IN))
    mexErrMsgTxt("Parameter `L' should be a structure.");
  if( (L_FIELD = mxGetField(L_IN,0,"perm")) == NULL)        /* L.perm */
    mexErrMsgTxt("Missing field L.perm.");
  if(m != mxGetM(L_FIELD) * mxGetN(L_FIELD))
    mexErrMsgTxt("L.perm size mismatches B");
  permPr = mxGetPr(L_FIELD);
  if( (L_FIELD = mxGetField(L_IN,0,"L")) == NULL)           /* L.L */
    mexErrMsgTxt("Missing field L.L.");
  if( m != mxGetM(L_FIELD) || m != mxGetN(L_FIELD) )
    mexErrMsgTxt("Size L.L mismatch.");
  if(!mxIsSparse(L_FIELD))
    mexErrMsgTxt("L.L should be sparse.");
  ljc = mxGetJc(L_FIELD);
  lir = mxGetIr(L_FIELD);
  if( (L_FIELD = mxGetField(L_IN,0,"xsuper")) == NULL)      /* L.xsuper */
    mexErrMsgTxt("Missing field L.xsuper.");
  nsuper = mxGetM(L_FIELD) * mxGetN(L_FIELD) - 1;
  if( nsuper > m )
    mexErrMsgTxt("Size L.xsuper mismatch.");
  xsuperPr = mxGetPr(L_FIELD);
/* ------------------------------------------------------------
   Allocate int-part of sparse output matrix X(m x n)
   Heuristically set nnz to nnz(B) + 4*m.
   ------------------------------------------------------------ */
  maxnnz = bjc[n] + 4 * m;
  xjc = (int *) mxCalloc(n + 1, sizeof(int));
  xir = (int *) mxCalloc(maxnnz, sizeof(int));
/* ------------------------------------------------------------
   Allocate working arrays:
   int invperm(m), snode(m), xsuper(nsuper+1), xlindx(nsuper+1), lindx(nnz(L)),
   iwork(nsuper).
   char cwork(nsuper).
   ------------------------------------------------------------ */
  invperm   = (int *) mxCalloc(m,sizeof(int)); 
  snode     = (int *) mxCalloc(m,sizeof(int)); 
  xsuper    = (int *) mxCalloc(nsuper+1,sizeof(int));
  xlindx    = (int *) mxCalloc(nsuper+1,sizeof(int));
  lindx     = (int *) mxCalloc(ljc[m], sizeof(int));
  iwork = (int *) mxCalloc(nsuper, sizeof(int));
  cwork = (char *) mxCalloc(nsuper, sizeof(char));
/* ------------------------------------------------------------
   Convert PERM, XSUPER to integer and C-Style
   ------------------------------------------------------------ */
  for(i = 0; i < m; i++){
    j = permPr[i];
    invperm[--j] = i;                /* so that invperm[perm[i]] = i */
  }
  for(i = 0; i <= nsuper; i++){
    j =  xsuperPr[i];
    xsuper[i] = --j;
  }
/* ------------------------------------------------------------
   Create "snode" from xsuper, and get compact subscript (xlindx,lindx)
   from (ljc,lir,xsuper), i.e. nz-pattern per supernode.
   ------------------------------------------------------------ */
  snodeCompress(xlindx,lindx,snode, ljc,lir,xsuper,nsuper);
/* ------------------------------------------------------------
   Compute nz structure after forward solve
   ------------------------------------------------------------ */
  symbfwmat(xjc, &xir, &maxnnz, bjc, bir, invperm, snode, xsuper,
            xlindx, lindx,
            nsuper, m, n, iwork, cwork);
/* ------------------------------------------------------------
   Create output matrix x
   ------------------------------------------------------------ */
  X_OUT = mxCreateSparse(m,n, 1,mxREAL);
  mxFree(mxGetJc(X_OUT));                    /* jc */
  mxFree(mxGetIr(X_OUT));                    /* ir */
  mxFree(mxGetPr(X_OUT));                    /* pr */
  xpr = (double *) mxCalloc(maxnnz,sizeof(double));
  mxSetJc(X_OUT, xjc);
  mxSetIr(X_OUT, xir);
  mxSetPr(X_OUT, xpr);
  mxSetNzmax(X_OUT, maxnnz);
  for(i = 0; i < maxnnz; i++)
    xpr[i] = 1.0;
/* ------------------------------------------------------------
   Release working arrays.
   ------------------------------------------------------------ */
  mxFree(cwork);
  mxFree(iwork);
  mxFree(lindx);
  mxFree(xlindx);
  mxFree(xsuper);
  mxFree(snode);
  mxFree(invperm);
}
