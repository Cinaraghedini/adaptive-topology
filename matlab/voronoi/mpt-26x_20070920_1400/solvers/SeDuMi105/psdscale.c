/*
%                                           y = psdscale(ud,x,K [,transp])
% PSDSCALE  Computes length lenud (=sum(K.s.^2)) vector y.
%   !transp (default) then y[k] = vec(Ldk' * Xk * Ldk)
%   transp == 1 then y[k] = vec(Udk' * Xk * Udk)
%   Uses pivot ordering ud.perm if available and nonempty.
%
%   SEE ALSO scaleK, factorK.
% **********  INTERNAL FUNCTION OF SEDUMI **********
function y = psdscale(ud,x,K [,transp])

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
#include "triuaux.h"
#include "blksdp.h"

/*    y = psdscale(ud,x,K [,transp]) */

#define Y_OUT plhs[0]
#define NPAROUT 1

#define UD_IN prhs[0]
#define X_IN prhs[1]
#define K_IN prhs[2]
#define NPARINMIN 3
#define TRANSP_IN prhs[3]
#define NPARIN 4

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
  const mxArray *UD_FIELD;
  int lenfull, lenud, sdplen, transp, fwsiz, i,k;
  double *fwork, *y,*permPr;
  const double *x,*ud;
  int *perm, *iwork;
  coneK cK;
  char use_pivot;

/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  if(nrhs < NPARIN){
    transp = 0;
    if(nrhs < NPARINMIN)
      mexErrMsgTxt("psdscale requires more input arguments.");
  }
  else
    transp = mxGetScalar(TRANSP_IN);
  if (nlhs > NPAROUT)
    mexErrMsgTxt("psdscale generates less output arguments.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Get statistics of cone K structure
   ------------------------------------------------------------ */
  lenud = cK.rDim + cK.hDim;
  lenfull = cK.lpN +  cK.qDim + lenud;
  sdplen = cK.rLen + cK.hLen;
/* ------------------------------------------------------------
   Get scale data: ud.{u,perm}.
   ------------------------------------------------------------ */
  if(!mxIsStruct(UD_IN)){
    if(mxGetM(UD_IN) * mxGetN(UD_IN) != lenud)          /* ud is vector */
      mexErrMsgTxt("ud size mismatch.");
    ud = mxGetPr(UD_IN);
    use_pivot = 0;
  }
  else{                         /* ud is structure */
    if( (UD_FIELD = mxGetField(UD_IN,0,"u")) == NULL)            /* ud.u */
      mexErrMsgTxt("Field ud.u missing.");
    if(mxGetM(UD_FIELD) * mxGetN(UD_FIELD) != lenud)
      mexErrMsgTxt("ud.u size mismatch.");
    ud = mxGetPr(UD_FIELD);
    UD_FIELD = mxGetField(UD_IN,0,"perm");                       /* ud.perm */
    if((use_pivot = (UD_FIELD != NULL))){
      if(mxGetM(UD_FIELD) * mxGetN(UD_FIELD) == sdplen)
        permPr = mxGetPr(UD_FIELD);
      else if(mxGetM(UD_FIELD) * mxGetN(UD_FIELD) == 0)
        use_pivot = 0;
      else
        mexErrMsgTxt("ud.perm size mismatch");
    }
  }
/* ------------------------------------------------------------
   Get input x
   ------------------------------------------------------------ */
  if(mxIsSparse(X_IN))
    mexErrMsgTxt("Sparse x not supported by this version of psdscale.");
  x = mxGetPr(X_IN);
/* ------------------------------------------------------------
   Validate x-input, and let x point to PSD part.
   ------------------------------------------------------------ */
  if(mxGetM(X_IN) * mxGetN(X_IN) == lenfull)
    x += cK.lpN + cK.qDim;
  else if(mxGetM(X_IN) * mxGetN(X_IN) != lenud)
    mexErrMsgTxt("size x mismatch.");
/* ------------------------------------------------------------
   Allocate output Y(lenud)
   ------------------------------------------------------------ */
  Y_OUT =  mxCreateDoubleMatrix(lenud, 1, mxREAL);
  y = mxGetPr(Y_OUT);
/* ------------------------------------------------------------
   Allocate fwork 2 * [ max(cK.rMaxn^2, 2*cK.hMaxn^2) ]
   iwork = int(sdplen)
   ------------------------------------------------------------ */
  fwsiz = MAX(SQR(cK.rMaxn),2*SQR(cK.hMaxn));
  fwork = (double *) mxCalloc( MAX(1,2 * fwsiz), sizeof(double));
  iwork = (int *) mxCalloc( MAX(1, sdplen), sizeof(int) );
/* ------------------------------------------------------------
   Convert Fortran to C-style in perm:
   ------------------------------------------------------------ */
  if(use_pivot){
    perm = iwork;
    for(k = 0; k < sdplen; k++){
      i = permPr[k];
      perm[k] = --i;
    }
  }
  else
    perm = (int *) NULL;
/* ------------------------------------------------------------
   The actual job is done here:.
   ------------------------------------------------------------ */
  psdscaleK(y, ud, perm, x, cK, transp, fwork);
/* ------------------------------------------------------------
   Release working arrays.
   ------------------------------------------------------------ */
  mxFree(fwork);
  mxFree(iwork);
}
