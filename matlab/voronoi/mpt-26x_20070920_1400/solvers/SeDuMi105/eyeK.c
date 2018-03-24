/*
    x = eyeK(K)
    yields identity solution w.r.t. symmetric cone K.

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

#include "mex.h"
#include "blksdp.h"

#define X_OUT plhs[0]
#define K_IN prhs[0]


/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   x = eyeK(K)
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
 int i,j,k, nk, lenfull;
 double *x;
 coneK cK;

/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
 if(nrhs < 1)
   mexErrMsgTxt("eyeK requires 1 input argument.");
 if(nlhs > 1)
   mexErrMsgTxt("eyeK generates 1 output argument.");
 /* ------------------------------------------------------------
    Disassemble cone K structure
    ------------------------------------------------------------ */
 conepars(K_IN, &cK);
 /* ------------------------------------------------------------
    Get statistics of cone K structure
    ------------------------------------------------------------ */
 lenfull = cK.lpN +  cK.qDim + cK.rDim + cK.hDim;
 for(k = 0; k < cK.rconeN; k++)
   lenfull += cK.rconeNL[k];
 /* ------------------------------------------------------------
    Allocate output x
    ------------------------------------------------------------ */
 X_OUT =  mxCreateDoubleMatrix(lenfull, 1, mxREAL);
 x = mxGetPr(X_OUT);
 /* ------------------------------------------------------------
    The actual job is done here:.
    ------------------------------------------------------------ */
 /* ------------------------------------------------------------
    LP: x = ones(K.l,1)
    ------------------------------------------------------------ */
 for(k = 0; k < cK.lpN; k++)
   x[k] = 1.0;
 x += cK.lpN;
 /* ------------------------------------------------------------
    LORENTZ: x(1) = sqrt(2)
    ------------------------------------------------------------ */
 for(k = 0; k < cK.lorN; k++){
   nk = cK.lorNL[k];
   x[0] = M_SQRT2;
   x += nk;
 }
 /* ------------------------------------------------------------
    RCONE: x(1) = x(2) = 1
    ------------------------------------------------------------ */
 for(k = 0; k < cK.rconeN; k++){
   nk = cK.rconeNL[k];
   x[0] = 1.0;
   x[1] = 1.0;
   x += nk;
 }
 /* ------------------------------------------------------------
    PSD:  x = eye(nk)
    ------------------------------------------------------------ */
 for(k = 0; k < cK.sdpN; k++){
   nk = cK.sdpNL[k];
   for(i = 0, j = 0; i < nk; i++, j += nk+1)
     x[j] = 1.0;
   x += (1 + (k >= cK.rsdpN)) * SQR(nk);
 }
}
