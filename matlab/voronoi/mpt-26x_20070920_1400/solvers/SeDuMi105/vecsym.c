/*
%                                                              y = vecsym(x,K)
% VECSYM  For the PSD submatrices, we let Yk = (Xk+Xk')/2
%
% SEE ALSO sedumi
% **********  INTERNAL FUNCTION OF SEDUMI **********
function y = vecsym(x,K)

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
#include "mex.h"
#include "blksdp.h"

#define Y_OUT plhs[0]

#define X_IN prhs[0]
#define K_IN prhs[1]

/* ************************************************************
   PROCEDURE symproj -- Y = (X+X')/2
   INPUT x, n - full n x n matrix x.
   OUTPUT y - on output, contains (x+x')/2
   ************************************************************ */
void symproj(double *y, const double *x, const int n)
{
  int colp,i,j;
  double yij;

  /* ------------------------------------------------------------
     x points to x(:,i);     x+colp = x(:,j).
     ------------------------------------------------------------ */
  for(i = 0; i < n; x += n, y += n, i++){
    y[i] = x[i];                         /* diagonal entry */
    for(colp = n + i, j=i+1; j<n; j++, colp += n){
      yij = (x[j] + x[colp]) / 2;         /* x(i,j)+x(j,i) */
      y[j] = yij;
      y[colp] = yij;
    }
  }
}

/* ************************************************************
   PROCEDURE skewproj -- Y = (X-X')/2
   INPUT x, n - full n x n matrix x.
   OUTPUT y - on output, contains (x-x')/2
   ************************************************************ */
void skewproj(double *y, const double *x, const int n)
{
  int colp,i,j;
  double yij;

  /* ------------------------------------------------------------
     x points to x(:,i);     x+colp = x(:,j).
     ------------------------------------------------------------ */
  for(i = 0; i < n; x += n, y += n, i++){
    y[i] = 0.0;                         /* diagonal entry */
    for(colp = n + i, j=i+1; j<n; j++, colp += n){
      yij = (x[j] - x[colp]) / 2;         /* x(j,i) - x(i,j) */
      y[j] = yij;
      y[colp] = -yij;                   /* conjugate */
    }
  }
}

/* ************************************************************
   PROCEDURE vecsymPSD - Let y = (x+x')/2.
   INPUT
     x     - length sum(K.s.^2) vector, to be symmetrized.
     rsdpN - number of real PSD blocks
     sdpN  - total number of PSD blocks, length(K.s)
     sdpNL - K.s, length sdpN array listing the orders of the PSD blocks.
   OUTPUT
     y - length sum(K.s.^2) vector, y = vecsym(x,K).
   ************************************************************ */
void vecsymPSD(double *y, const double *x,const int rsdpN,const int sdpN,
               const double *sdpNL)
{
  int k,nk,nksqr;
/* ------------------------------------------------------------
   Make real PSD blocks symmetric
   ------------------------------------------------------------ */
  for(k = 0; k < rsdpN; k++){
    nk = sdpNL[k]; nksqr = SQR(nk);
    symproj(y,x,nk);
    x += nksqr; y += nksqr;
  }
/* ------------------------------------------------------------
   Make complex PSD blocks Hermitian
   ------------------------------------------------------------ */
  for(; k < sdpN; k++){
    nk = sdpNL[k]; nksqr = SQR(nk);
    symproj(y,x,nk);
    x += nksqr; y += nksqr;
    skewproj(y,x,nk);
    x += nksqr; y += nksqr;
  }
}

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
     y = vecsym(x,K)

     Computes "symmetrization of x: Yk = (Xk+Xk')/2
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
 mxArray *output_array[1], *Xk;

 coneK cK;
 int k, nk, nksqr, lqDim,lenfull;
 const double *x;
 double *y;

 /* ------------------------------------------------------------
    Check for proper number of arguments
    ------------------------------------------------------------ */
 if(nrhs < 2)
   mexErrMsgTxt("vecsym requires 2 input arguments.");
 if(nlhs > 1)
   mexErrMsgTxt("vecsym generates 1 output argument.");
 /* ------------------------------------------------------------
    Disassemble cone K structure
    ------------------------------------------------------------ */
 conepars(K_IN, &cK);
 /* ------------------------------------------------------------
    Compute some statistics based on cone K structure
    ------------------------------------------------------------ */
 lqDim = cK.lpN + cK.qDim;
 lenfull = lqDim + cK.rDim + cK.hDim;
 /* ------------------------------------------------------------
    Get input vector x
    ------------------------------------------------------------ */
 if(mxIsSparse(X_IN))
   mexErrMsgTxt("x must be full.");
 if(mxGetM(X_IN) * mxGetN(X_IN) != lenfull) 
   mexErrMsgTxt("Parameter `x' size mismatch.");
 x = mxGetPr(X_IN);
/* ------------------------------------------------------------
   Allocate output vector y, and make it vecsym(x)
   ------------------------------------------------------------ */
 Y_OUT = mxCreateDoubleMatrix(lenfull, 1, mxREAL);
 y = mxGetPr(Y_OUT);
 memcpy(y,x,lqDim * sizeof(double));
 x += lqDim; y += lqDim;
 vecsymPSD(y, x,cK.rsdpN,cK.sdpN,cK.sdpNL);
}
