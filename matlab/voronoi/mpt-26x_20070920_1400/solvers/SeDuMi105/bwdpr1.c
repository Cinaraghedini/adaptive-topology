/*
%                                                            y = bwdpr1(Lden, b)
% BWDPR1  Solves "PROD_k L(pk,betak)' * y = b", where
%     L(p,beta) = eye(n) + tril(p*beta',-1).
%
% SEE ALSO sedumi, dpr1fact, fwdpr1
% **********  INTERNAL FUNCTION OF SEDUMI **********
function y = bwdpr1(Lden, b)

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

#include <math.h>
#include <string.h>
#include "mex.h"
#include "blksdp.h"

/* y = bwdpr1fact(Lden,b) */
#define Y_OUT plhs[0]
#define NPAROUT 1

#define LDEN_IN prhs[0]
#define B_IN prhs[1]
#define NPARIN 2

/* ------------------------------------------------------------
   PROCEDURE bwipr1 - I(dentity) P(lus) R(ank)1 forward solve.
   INPUT:
   p    - length m floats
   beta - length n floats
   m, n - order of p and beta, resp. (n <= m)
   UPDATED:
   y - Length m. On input, contains the rhs. On output, the solution to
       L(p,beta)'*yNEW = yOLD. This updates only y(0:n-1).
   ------------------------------------------------------------ */
void bwipr1(double *y, const double *p, const double *beta,
            const int m, const int n)
{
  int i;
  double yi,t;

  if(n < 1)           /* If L = I, y remains the same */
    return;
/* ------------------------------------------------------------
   Let t = p(n:m-1)' * y
   ------------------------------------------------------------ */
  t = realdot(p+n, y+n, m-n);
/* ------------------------------------------------------------
   Solve yi for i = n-1:-1:0, from
   (eye(m)+triu(beta*p',1)) * yNEW = yOLD,
   i.e. yNEW(i) + betai*t = yOLD(i), with t := p(i+1:n-1)'*y.
   ------------------------------------------------------------ */
  for(i = n-1; i >= 0; i--){
    yi = (y[i] -= t * beta[i]);
    t += p[i] * yi;
  }
}

/* ------------------------------------------------------------
   PROCEDURE bwipr1o - I(dentity) P(lus) R(ank)1 forward solve, O(rdered).
   INPUT:
   perm - length m permutation on p and y.
   p    - length m floats
   beta - length n floats
   m, n - order of p and beta, resp. (n <= m)
   UPDATED:
   y - Length m. On input, contains the rhs. On output, the solution to
       L(p,beta)'*yNEW = yOLD. This updates only y(0:n-1).
   ------------------------------------------------------------ */
void bwipr1o(double *y, const int *perm, const double *p, const double *beta,
            const int m, const int n)
{
  int i, permi;
  double yi,t;

  if(n < 1)           /* If L = I, y remains the same */
    return;
/* ------------------------------------------------------------
   Let t = p(perm(n:m-1))' * y(perm(n:m-1))
   ------------------------------------------------------------ */
  for(t = 0.0, i = m-1; i >= n; i--){
    permi = perm[i];
    t += p[permi] * y[permi];
  }
/* ------------------------------------------------------------
   Solve yi for i = n-1:-1:0, from
   (eye(m)+triu(beta*p',1)) * yNEW = yOLD,
   i.e. yNEW(i) + betai*t = yOLD(i), with t := p(i+1:n-1)'*y.
   ------------------------------------------------------------ */
  for(; i >= 0; i--){
    permi = perm[i];
    yi = (y[permi] -= t * beta[i]);
    t += p[permi] * yi;
  }
}

/* ************************************************************
   PROCEDURE bwprodform - Solves (PROD_j L(pj,betaj))' * yNEW = yOLD.
   INPUT
     p - nonzeros of sparse m x n matrix P. Has xsuper(j+1) nonzeros in
      column j.
     xsuper - xsuper(j+1) is number of nonzeros in p(:,j).
     perm - lists pivot order for columns where ordered(j)==1.
     ordered - ordered[j]==1 iff p(:,j) and beta(L:,j) have been reordered;
       the original row numbers are in perm(:,j).
     n - number of columns in p, beta.
   UPDATED
     y - length m vector. On input, the rhs. On output the solution to
       (PROD_j L(pj,betaj))' * yNEW = yOLD.
   ************************************************************ */
void bwprodform(double *y, const int *xsuper, const int *perm,
                const double *p, const double *beta, const int *betajc,
                const char *ordered, const int n, int pnnz,
                int permnnz)
{
  int k,nk, mk;
/* ------------------------------------------------------------
   Backward solve L(pk,betak) * yNEXT = yPREV   for k=n-1:-1:0.
   ------------------------------------------------------------ */
  for(k = n-1; k >= 0; k--){
    mk = xsuper[k+1];
    nk = betajc[k+1] - betajc[k];
    pnnz -= mk;
    if(ordered[k]){
      permnnz -= mk;
      bwipr1o(y, perm+permnnz, p+pnnz, beta+betajc[k], mk, nk);
    }
    else
      bwipr1(y, p+pnnz, beta+betajc[k], mk, nk);
  }
  mxAssert(pnnz == 0,"");
  mxAssert(permnnz == 0 || permnnz == 1,"");
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
 char *ordered;
 int m,n,nden,dznnz, i,j, permnnz, pnnz;
 const double *beta, *betajcPr, *orderedPr, *pivpermPr, *p;
 int *betajc, *pivperm;
 double *y, *fwork;
 jcir dz;
/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  if(nrhs < NPARIN)
    mexErrMsgTxt("bwdpr1 requires more input arguments.");
  if(nlhs > NPAROUT)
    mexErrMsgTxt("bwdpr1 generates less output arguments.");
/* ------------------------------------------------------------
   Get input b
   ------------------------------------------------------------ */
  if(mxIsSparse(B_IN))
    mexErrMsgTxt("b should be full");
  m = mxGetM(B_IN);
  n = mxGetN(B_IN);
/* ------------------------------------------------------------
   Create output y as a copy of b
   ------------------------------------------------------------ */
  Y_OUT = mxDuplicateArray(B_IN);                 /* Y_OUT = B_IN */
  y = mxGetPr(Y_OUT);
/* ------------------------------------------------------------
   Disassemble dense-update structure Lden (p,xsuper,beta,betajc,rowperm)
   NOTE: if there are no dense columns, then simply let y = b, and return.
   ------------------------------------------------------------ */
  if(!mxIsStruct(LDEN_IN))
    mexErrMsgTxt("Parameter `Lden' should be a structure.");
  if( (MY_FIELD = mxGetField(LDEN_IN,0,"betajc")) == NULL)      /* betajc */
    mexErrMsgTxt("Missing field Lden.betajc.");
  nden = mxGetM(MY_FIELD) * mxGetN(MY_FIELD) - 1;
/* If no dense columns: get out immediately ! */
  if(nden == 0)
    return;
  else{
    betajcPr = mxGetPr(MY_FIELD);
    if( (MY_FIELD = mxGetField(LDEN_IN,0,"p")) == NULL)         /* p */
      mexErrMsgTxt("Missing field Lden.p.");
    p = mxGetPr(MY_FIELD);
    pnnz = mxGetM(MY_FIELD) * mxGetN(MY_FIELD);
    if( (MY_FIELD = mxGetField(LDEN_IN,0,"dopiv")) == NULL)     /* dopiv */
      mexErrMsgTxt("Missing field Lden.dopiv.");
    if(mxGetM(MY_FIELD) * mxGetN(MY_FIELD) != nden)
      mexErrMsgTxt("Size mismatch Lden.dopiv.");
    orderedPr = mxGetPr(MY_FIELD);
    if( (MY_FIELD = mxGetField(LDEN_IN,0,"pivperm")) == NULL)   /* pivperm */
      mexErrMsgTxt("Missing field Lden.pivperm.");
    pivpermPr = mxGetPr(MY_FIELD);
    permnnz = mxGetM(MY_FIELD) * mxGetN(MY_FIELD);
    if( (MY_FIELD = mxGetField(LDEN_IN,0,"beta")) == NULL)      /* beta */
      mexErrMsgTxt("Missing field Lden.beta.");
    beta = mxGetPr(MY_FIELD);
    if( (MY_FIELD = mxGetField(LDEN_IN,0,"dz")) == NULL)         /* dz */
      mexErrMsgTxt("Missing field Lden.dz.");
    if(mxGetM(MY_FIELD) != m || mxGetN(MY_FIELD) != nden)
      mexErrMsgTxt("Lden.dz size mismatch.");
    if(!mxIsSparse(MY_FIELD))
      mexErrMsgTxt("Lden.dz must be sparse.");
    dz.jc = mxGetJc(MY_FIELD);
    dz.ir = mxGetIr(MY_FIELD);                             /* (rowperm) */
    dznnz = dz.jc[nden];
    if(dznnz > m)
      mexErrMsgTxt("Lden.dz size mismatch.");
/* ------------------------------------------------------------
   Allocate working arrays int: betajc(nden+1), pivperm(permnnz),
   char: ordered(nden)
   double: fwork(dznnz)
   ------------------------------------------------------------ */
    betajc = (int *) mxCalloc(nden + 1,sizeof(int));
    pivperm = (int *) mxCalloc(MAX(permnnz,1),sizeof(int));
    ordered = (char *) mxCalloc(nden,sizeof(char));   /* nden > 0 */
    fwork = (double *) mxCalloc(MAX(dznnz,1), sizeof(double));
/* ------------------------------------------------------------
   Convert betajcPr, ordered, pivperm to int
   ------------------------------------------------------------ */
    for(i = 0; i <= nden; i++){
      j = betajcPr[i];
      betajc[i] = --j;
    }
    for(i = 0; i < nden; i++){
      ordered[i] = orderedPr[i];
    }
    for(i = 0; i < permnnz; i++){
      pivperm[i] = pivpermPr[i];
    }
    MY_FIELD = mxGetField(LDEN_IN,0,"beta");
    if(mxGetM(MY_FIELD) * mxGetN(MY_FIELD) != betajc[nden])
      mexErrMsgTxt("Size mismatch Lden.beta.");
/* ------------------------------------------------------------
   The actual job is done here: y(perm,:) = PROD_L'\b(perm,:).
   ------------------------------------------------------------ */
    for(j = 0; j < n; j++, y += m){
      for(i = 0; i < dznnz; i++)            /* fwork = y(dzir) */
        fwork[i] = y[dz.ir[i]];
      bwprodform(fwork, dz.jc, pivperm, p, beta, betajc, ordered, nden,
                 pnnz, permnnz);
      for(i = 0; i < dznnz; i++)            /* y(dzir) = fwork */
        y[dz.ir[i]] = fwork[i];
    }
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
    mxFree(fwork);
    mxFree(ordered);
    mxFree(pivperm);
    mxFree(betajc);
  } /* if dense columns */
}
