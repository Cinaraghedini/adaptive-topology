/*
%                                                y = qblkmul(mu,d,blkstart)
% QBLKMUL  yields length(y)=blkstart(end)-blkstart(1) vector with
%    y[k] = mu(k) * d[k]; the blocks d[k] are partitioned by blkstart.
%
% SEE ALSO sedumi
% **********  INTERNAL FUNCTION OF SEDUMI **********
function y = qblkmul(mu,d,blkstart)

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

#include "mex.h"
#include "blksdp.h"

#define Y_OUT plhs[0]
#define NPAROUT 1
#define MU_IN prhs[0]
#define D_IN prhs[1]
#define BLKSTART_IN prhs[2]
#define NPARIN 3

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
  int i,j,nblk,k,nk,qDim;
  double *y;
  const double *d, *mu, *blkstartPr;
  int *blkstart;
/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  if(nrhs < NPARIN)
    mexErrMsgTxt("qblkmul requires more input arguments.");
  if (nlhs > NPAROUT)
    mexErrMsgTxt("qblkmul generates 1 output argument.");
/* ------------------------------------------------------------
   Get inputs d, mu, blkstart
   ------------------------------------------------------------ */
  d = mxGetPr(D_IN);
  qDim = mxGetM(D_IN) * mxGetN(D_IN);
  nblk = mxGetM(MU_IN) * mxGetN(MU_IN);
  mu = mxGetPr(MU_IN);
  blkstartPr = mxGetPr(BLKSTART_IN);
  if(nblk != mxGetM(BLKSTART_IN) * mxGetN(BLKSTART_IN) - 1)
    mexErrMsgTxt("blkstart size mismatch.");
/* ------------------------------------------------------------
   Allocate int working array blkstart(nblk+1).
   ------------------------------------------------------------ */
  blkstart = (int *) mxCalloc(nblk + 1, sizeof(int));
/* ------------------------------------------------------------
   Convert Fortran double to C int
   ------------------------------------------------------------ */
  for(i = 0; i <= nblk; i++){
    j = blkstartPr[i];             /* double to int */
    blkstart[i] = --j;
  }
/* ------------------------------------------------------------
   Let d point to Lorentz norm-bound
   ------------------------------------------------------------ */
  if(qDim != blkstart[nblk] - blkstart[0]){
    if(qDim == nblk + blkstart[nblk] - blkstart[0]){
      d += nblk;                   /* Point to Lorentz norm-bound */
    }
    else if(qDim >= blkstart[nblk]){
      d += blkstart[0];                   /* Point to Lorentz norm-bound */
    }
    else
      mexErrMsgTxt("d size mismatch.");
  }
  qDim = blkstart[nblk] - blkstart[0];
/* ------------------------------------------------------------
   Allocate output y(qDim)
   ------------------------------------------------------------ */
  Y_OUT =  mxCreateDoubleMatrix(qDim, 1, mxREAL);
  y = mxGetPr(Y_OUT);
/* ------------------------------------------------------------
   LORENTZ: yk = mu(k) * d[k]
   ------------------------------------------------------------ */
  for(k = 0; k < nblk; k++){
    nk = blkstart[k+1] - blkstart[k];
    scalarmul(y,mu[k],d,nk);
    y += nk;
    d += nk;
  }
}
