/* ************************************************************
   MODULE sdmaux*.c  -- Several low-level subroutines for the
   mex-files in the Self-Dual-Minimization package.

    This file is part of SeDuMi 1.04BETA
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
#include "blksdp.h"
/* ============================================================
   CONE-K STATISTICS
   ============================================================ */
/* ************************************************************
   PROCEDURE conepars - Read cone K parameters from K-structure
   INPUT
     mxK  -  the Matlab structure "K", as passes as input argument "K_IN".
   OUTPUT
     *pK - struct where cone K parameters get stored.
   ************************************************************ */
void conepars(const mxArray *mxK, coneK *pK)
{
 const mxArray *K_FIELD;
 const double *blkstartPr;
 int idummy, nblk;
 char gotthem;

 if(!mxIsStruct(mxK))
   mexErrMsgTxt("Parameter `K' should be a structure.");
 if( (K_FIELD = mxGetField(mxK,0,"f")) == NULL)      /* K.f */
   pK->frN = 0;
 else
   pK->frN = mxGetScalar(K_FIELD);
 if( (K_FIELD = mxGetField(mxK,0,"l")) == NULL)      /* K.l */
   pK->lpN = 0;
 else
   pK->lpN = mxGetScalar(K_FIELD);
 if( (K_FIELD = mxGetField(mxK,0,"q")) == NULL)      /* K.q */
   pK->lorN = 0;
 else{
   pK->lorN = mxGetM(K_FIELD) * mxGetN(K_FIELD);
   pK->lorNL = mxGetPr(K_FIELD);
   if(pK->lorN == 1)                                /* K.q=0 -> lorN = 0 */
     if(pK->lorNL[0] == 0.0)
       pK->lorN = 0;
 }
 if( (K_FIELD = mxGetField(mxK,0,"r")) == NULL)      /* K.r */
   pK->rconeN = 0;
 else{
   pK->rconeN = mxGetM(K_FIELD) * mxGetN(K_FIELD);
   pK->rconeNL = mxGetPr(K_FIELD);
   if(pK->rconeN == 1)                                /* K.r=0 -> rconeN = 0 */
     if(pK->rconeNL[0] == 0.0)
       pK->rconeN = 0;
 }
 if( (K_FIELD = mxGetField(mxK,0,"s")) == NULL){     /* K.s */
   pK->sdpN = 0;
 }
 else{
   pK->sdpN = mxGetM(K_FIELD) * mxGetN(K_FIELD);
   pK->sdpNL = mxGetPr(K_FIELD);
   if(pK->sdpN == 1)                                /* K.s=0 -> sdpN = 0 */
     if(pK->sdpNL[0] == 0.0)
       pK->sdpN = 0;
 }
 if( (K_FIELD = mxGetField(mxK,0,"rsdpN")) == NULL)      /* K.rsdpN */
   pK->rsdpN = pK->sdpN;                           /* default to all real */
 else
   if((pK->rsdpN = mxGetScalar(K_FIELD)) > pK->sdpN)
     mexErrMsgTxt("K.rsdpN mismatches K.s");
 /* --------------------------------------------------
    GET STATISTICS: try to read from K, otherwise compute them.
    -------------------------------------------------- */
 gotthem = 0;
 if( (K_FIELD = mxGetField(mxK,0,"rLen")) != NULL){      /* K.rLen */
   pK->rLen = mxGetScalar(K_FIELD);
   if( (K_FIELD = mxGetField(mxK,0,"hLen")) != NULL){      /* K.hLen */
     pK->hLen = mxGetScalar(K_FIELD);
     if( (K_FIELD = mxGetField(mxK,0,"qMaxn")) != NULL){      /* K.qMaxn */
       pK->qMaxn = mxGetScalar(K_FIELD);
       if( (K_FIELD = mxGetField(mxK,0,"rMaxn")) != NULL){      /* K.rMaxn */
	 pK->rMaxn = mxGetScalar(K_FIELD);
	 if( (K_FIELD = mxGetField(mxK,0,"hMaxn")) != NULL){    /* K.hMaxn */
	   pK->hMaxn = mxGetScalar(K_FIELD);
	   if( (K_FIELD = mxGetField(mxK,0,"blkstart"))!=NULL){ /*K.blkstart*/
	     if(mxIsSparse(K_FIELD))
	       mexErrMsgTxt("K.blkstart must be a full vector.");
             nblk = 1 + pK->lorN + pK->sdpN;
             if(mxGetM(K_FIELD) * mxGetN(K_FIELD) != nblk + 1)
	       mexErrMsgTxt("Size mismatch K.blkstart.");
	     blkstartPr = mxGetPr(K_FIELD);
	     pK->qDim = blkstartPr[pK->lorN+1] - blkstartPr[0];
	     blkstartPr += pK->lorN+1;
	     pK->rDim = blkstartPr[pK->rsdpN] - blkstartPr[0];
	     pK->hDim = blkstartPr[pK->sdpN] - blkstartPr[pK->rsdpN];
	     gotthem = 1;
	   }
	 }
       }
     }
   }
 }
 if(!gotthem){
   someStats(&(pK->qMaxn), &(pK->qDim), &idummy, pK->lorNL, pK->lorN);
   someStats(&(pK->rMaxn), &(pK->rLen), &(pK->rDim), pK->sdpNL, pK->rsdpN);
   someStats(&(pK->hMaxn), &(pK->hLen), &(pK->hDim), pK->sdpNL+pK->rsdpN,
	     (pK->sdpN) - (pK->rsdpN));
   pK->hDim *= 2;
 }
}

/* ************************************************************
   PROCEDURE someStats  --  Computes maximum, sum and sum of squares
   INPUT
   x, n - length n vector
   OUTPUT
   xmax, xsum, xssqr - Maximum, sum total and sum of squares
   IMPORTANT: this routine is especially designed for use with the
    blk.s structure, which contains nonneg integers stored as doubles.
   ************************************************************ */
void someStats(int *pxmax, int *pxsum, int *pxssqr,
	       const double *x, const int n)
{
 int xi, xmax, xsum, xssqr;
 int i;

 xmax = 0;             /* assume that all integers are nonnegative */
 xsum = 0; xssqr = 0;
 for(i = 0; i < n; i++){
   xi = x[i];
   xmax = MAX(xmax, xi);
   xsum += xi;
   xssqr += SQR(xi);
 }
 *pxmax = xmax;
 *pxsum = xsum;
 *pxssqr = xssqr;
}
