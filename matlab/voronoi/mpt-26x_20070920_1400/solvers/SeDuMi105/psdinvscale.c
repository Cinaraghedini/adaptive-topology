/*
   y = psdinvscale(ud,x,K,transp)
   Computes y = D(d^{-1}) x with d in K.
   !transp then Y = Ld' \ X / Ld
   transp == 1 then Y = Ud' \ X / Ud

    This file is part of SeDuMi 1.03BETA
    Copyright (C) 1999 Jos F. Sturm
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
#include "triuaux.h"
#include "blksdp.h"

/*    y = psdinvscale(ud,x,K,transp) */

#define Y_OUT plhs[0]
#define NPAROUT 1

#define UD_IN prhs[0]
#define X_IN prhs[1]
#define K_IN prhs[2]
#define NPARINMIN 3
#define TRANSP_IN prhs[3]
#define NPARIN 4

/* ************************************************************
   PROCEDURE psdinvscale - Computes y = D(d^{-1})x.
   (transp == 0) Y = U\X/U'
   (transp == 1) Y = U'\X/U.
   INPUT
     x - length lenud input vector.
     ud - Cholesky factor of d for PSD part (after PERM ordering).
     psdNL - length psdN array
     rpsdN - number of real symmetric PSD blocks
     sdpN  - total number of psd blocks
     transp - boolean
   OUTPUT
     y - length lenud output vector, y=D(d)x.
   WORK
     fwork - length max(rmaxn^2,2*hmaxn^2) working vector.
   ************************************************************ */
void psdinvscale(double *y, const double *ud, const double *x,
                 const int *psdNL, const int rsdpN, const int sdpN,
                 const char transp, double *fwork)
{
  int k,nk,nksqr;
/* ------------------------------------------------------------
   PSD, !transp: triu(Y) = triu(Ld' \ (X(perm,perm) / Ld)).
   Needs ony tril(X/Ld).
   ------------------------------------------------------------ */
  if(!transp){
    for(k = 0; k < rsdpN; k++){                /* real symmetric */
      nk = psdNL[k];
      invltxl(y,ud,x,nk,fwork);
      tril2sym(y,nk);
      nksqr = SQR(nk);
      y += nksqr; ud += nksqr;
      x += nksqr;
    }
    for(; k < sdpN; k++){                    /* complex Hermitian */
      nk = psdNL[k];
      nksqr = SQR(nk);
      prpiinvltxl(y,y+nksqr,ud,ud+nksqr,x,x+nksqr,nk,fwork);
      tril2herm(y,y+nksqr,nk);
      nksqr += nksqr;
      y += nksqr; ud += nksqr;
      x += nksqr;
    }
  } /* !transp */
  else{
/* ------------------------------------------------------------
   PSD, transp: triu(Y) = triu(Ud' \ (X(perm,perm) / Ud)).
   Needs ony triu(X/Ud).
   ------------------------------------------------------------ */
    for(k = 0; k < rsdpN; k++){                /* real symmetric */
      nk = psdNL[k];
      invutxu(y,ud,x,nk,fwork);
      triu2sym(y,nk);
      nksqr = SQR(nk);
      y += nksqr; ud += nksqr;
      x += nksqr;
    }
    for(; k < sdpN; k++){                    /* complex Hermitian */
      nk = psdNL[k];
      nksqr = SQR(nk);
      prpiinvutxu(y,y+nksqr,ud,ud+nksqr,x,x+nksqr,nk,fwork);
      triu2herm(y,y+nksqr,nk);
      nksqr += nksqr;
      y += nksqr; ud += nksqr;
      x += nksqr;
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
  int lenfull, lenud, transp, fwsiz, i, ifirst;
  double *fwork, *y;
  const double *x,*ud;
  int *psdNL;
  coneK cK;
/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  if(nrhs < NPARIN){
    transp = 0;
    if(nrhs < NPARINMIN)
      mexErrMsgTxt("psdinvscale requires more input arguments.");
  }
  else
    transp = mxGetScalar(TRANSP_IN);
  if (nlhs > NPAROUT)
    mexErrMsgTxt("psdinvscale generates 1 output argument.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Get statistics of cone K structure
   ------------------------------------------------------------ */
  ifirst = cK.lpN +  cK.qDim;
  lenud = cK.rDim + cK.hDim;
  lenfull = ifirst + lenud;
/* ------------------------------------------------------------
   Get inputs ud, x.
   ------------------------------------------------------------ */
  if(mxGetM(UD_IN) * mxGetN(UD_IN) != lenud)
    mexErrMsgTxt("ud size mismatch.");
  ud = mxGetPr(UD_IN);
  if(mxIsSparse(X_IN))
    mexErrMsgTxt("Sparse x not supported by this version of psdinvscale.");
  if(mxGetM(X_IN) * mxGetN(X_IN) != lenfull)
    mexErrMsgTxt("size x mismatch.");
  x = mxGetPr(X_IN) + ifirst;           /* Jump to PSD part */
/* ------------------------------------------------------------
   Allocate output Y
   ------------------------------------------------------------ */
  Y_OUT =  mxCreateDoubleMatrix(lenud, 1, mxREAL);
  y = mxGetPr(Y_OUT);
/* ------------------------------------------------------------
   Allocate fwork [ max(cK.rMaxn^2, 2*cK.hMaxn^2) ]
   psdNL = int(cK.sdpN)
   ------------------------------------------------------------ */
  fwsiz = MAX(SQR(cK.rMaxn),2*SQR(cK.hMaxn));
  fwork = (double *) mxCalloc( MAX(1,fwsiz), sizeof(double));
  psdNL = (int *) mxCalloc( MAX(1, cK.sdpN), sizeof(int) );
/* ------------------------------------------------------------
   Convert double to int
   ------------------------------------------------------------ */
  for(i = 0; i < cK.sdpN; i++)
    psdNL[i] = cK.sdpNL[i];               /* double to int */
/* ------------------------------------------------------------
   The real job:
   ------------------------------------------------------------ */
  psdinvscale(y, ud,x,psdNL,cK.rsdpN,cK.sdpN,transp, fwork);
/* ------------------------------------------------------------
   Release working arrays.
   ------------------------------------------------------------ */
  mxFree(fwork);
  mxFree(psdNL);
}
