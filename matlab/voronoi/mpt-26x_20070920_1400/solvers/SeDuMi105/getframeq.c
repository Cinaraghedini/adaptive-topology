/*
    frm = getframeq(x,lab, K)

    This file is part of SeDuMi 1.05 DEVELOPMENT
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

#include <math.h>       /* for def sqrt() */
#include "mex.h"
#include "blksdp.h"

#define FRM_OUT plhs[0]
#define NPAROUT 1
#define X_IN prhs[0]
#define LAB_IN prhs[1]
#define K_IN prhs[2]
#define NPARIN 3

/* ************************************************************
   PROCEDURE getframeq
   INPUT
     x - length qDim=sum(lorNL)
     lab - 2*lorN
     lorNL - lorN
     lorN - number of Lorentz blocks
    OUTPUT
      frm - size sum(lorNL)-lorN. kth block is length lorNL[k]-1 vector
        of norm 1/sqrt(2).
   ************************************************************ */
#define GETFRAME_TOL 1E-10
void getframeq(double *frm, const double *x, const double *lab,
               const int *lorNL, const int lorN)
{
  int k,nk;
  double difk;
  const double *lab2;
/* ------------------------------------------------------------
   Let lab2 = lab + lorN point to the large spectral value.
   ------------------------------------------------------------ */
 lab2 = lab + lorN;
 ++x;
/* ------------------------------------------------------------
   for each k: frmk = xk(2:nk) / (labk(2)-labk(1))
   ------------------------------------------------------------ */
 for(k = 0; k < lorN; k++){
   difk = lab2[k] - lab[k];
   nk = lorNL[k];
   if(difk > GETFRAME_TOL * lab2[k])
     scalarmul(frm, 1/difk, x, nk-1);
   else{
/* ------------------------------------------------------------
   If (lab2-lab1)/lab2 is almost zero, then compute ||xk(2:nk)||
   from scratch. If it is exactly zero, then take arbitrary frame.
   ------------------------------------------------------------ */
     difk = realssqr(x,nk-1);
     if(difk > 0)
       scalarmul(frm, 1/sqrt(2*difk), x, nk-1);  /* norm = 1/sqrt(2) */
     else{
       frm[0] = M_SQRT1_2;        /* identity -> take e1/sqrt2 frame */
       fzeros(frm+1, nk-2);
     }
   }
/* ------------------------------------------------------------
   Point to next block
   ------------------------------------------------------------ */
   x += nk;
   frm += nk-1;
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
 int i;
 double *frm;
 const double *x, *lab;
 coneK cK;
 int *iwork;

/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
 if(nrhs < NPARIN)
   mexErrMsgTxt("getframeq requires more input arguments.");
 if (nlhs > NPAROUT)
   mexErrMsgTxt("getframeq generates 1 output argument.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
 conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Get inputs x, lab
   ------------------------------------------------------------ */
 x = mxGetPr(X_IN);
 if(mxGetM(X_IN) * mxGetN(X_IN) != cK.qDim)
   if(mxGetM(X_IN) * mxGetN(X_IN) < cK.lpN + cK.qDim)
     mexErrMsgTxt("x size mismatch");
   else
     x += cK.lpN;                   /* Point to Lorentz part */
 lab = mxGetPr(LAB_IN);
 if(mxGetM(LAB_IN) * mxGetN(LAB_IN) != 2*cK.lorN)
   if(mxGetM(LAB_IN) * mxGetN(LAB_IN) < cK.lpN + 2*cK.lorN)
     mexErrMsgTxt("lab size mismatch");
   else
     lab += cK.lpN;
/* ------------------------------------------------------------
   Allocate output frm(qDim-lorN)
   ------------------------------------------------------------ */
 FRM_OUT =  mxCreateDoubleMatrix(cK.qDim-cK.lorN, 1, mxREAL);
 frm = mxGetPr(FRM_OUT);
/* ------------------------------------------------------------
   Allocate working array iwork(lorN), and let iwork = lorNL.
   ------------------------------------------------------------ */
 iwork = (int *) mxCalloc(MAX(cK.lorN,1), sizeof(int));
 for(i = 0; i < cK.lorN; i++)
   iwork[i] = cK.lorNL[i];           /* float to int */
/* ------------------------------------------------------------
   The real job: let frmk = xk(2:nk)/(lab2(k)-lab1(k)).
   ------------------------------------------------------------ */
 getframeq(frm, x, lab, iwork, cK.lorN);
/* ------------------------------------------------------------
   RELEASE WORKING ARRAY
   ------------------------------------------------------------ */
 mxFree(iwork);
}
