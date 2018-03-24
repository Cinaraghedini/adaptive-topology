/* ************************************************************
%                                                K = statsK(K)
% STATSK  Collects statistics (max and sum of dimensions) of cone K
%
% SEE ALSO sedumi.
% ******************** INTERNAL FUNCTION OF SEDUMI ********************
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

   ************************************************************ */

#include <string.h>
#include "mex.h"
#include "blksdp.h"

/* K = statsK(K) */

#define K_OUT plhs[0]
#define NPAROUT 1

#define K_IN prhs[0]
#define NPARIN 1

/* ************************************************************
   PROCEDURE setScalarField
   INPUT
     fnm  -
     fval -
   UPDATED
     mxK  -
   ************************************************************ */
void setScalarField(mxArray *mxK,const char *fnm,const double fval)
{
  mxArray *K_FIELD;
  K_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL);
  *mxGetPr(K_FIELD) = fval;
  mxSetField(mxK,0,fnm,K_FIELD);
}

/* ************************************************************
   PROCEDURE getblkstart - computes 2+|K.q|+|K.s| blkstart vector.
   INPUT
     cK - cone structure K.{l,q,s}.
     blkstart - blkstart = cumsum([0 K.l K.q K.s.^2])
   ************************************************************ */
void getblkstart(int *blkstart, const coneK cK)
{
  int i,k,nk,blknz;

/* ------------------------------------------------------------
   Lorentz trace block
   ------------------------------------------------------------ */
  blkstart[0] = cK.lpN;
  blknz = cK.lpN + cK.lorN;
  blkstart[1] = blknz;
  k = 1; 
/* ------------------------------------------------------------
   Lorentz norm-bound blocks
   ------------------------------------------------------------ */
  for(i = 0; i < cK.lorN; i++)
    blkstart[++k] = (blknz += cK.lorNL[i]-1);
/* ------------------------------------------------------------
   PSD block
   ------------------------------------------------------------ */
  for(i = 0; i < cK.rsdpN; i++){           /* real symmetric */
    nk = cK.sdpNL[i];
    blkstart[++k] = (blknz += SQR(nk));
  }
  for(; i < cK.sdpN; i++){                  /* complex Hermitian */
    nk = cK.sdpNL[i];
    blkstart[++k] = (blknz += 2*SQR(nk));
  }
}

#define NKFIELDS 11
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
  int nblk,i, nycomplex;
  int *blkstart;
  double *blkstartPr;
  const double *ycomplexPr;
  coneK cK;
  mxArray *K_FIELD;
  const char *Kfieldnames[] = {"l",    "q",    "s",     "rsdpN", "blkstart",
                               "rLen", "hLen", "qMaxn", "rMaxn", "hMaxn",
                               "ycomplex"};
/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  if(nrhs < NPARIN)
    mexErrMsgTxt("statsK requires more input arguments.");
  if (nlhs > NPAROUT)
    mexErrMsgTxt("statsK generates less output arguments.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
  if( (K_FIELD = mxGetField(K_IN,0,"ycomplex")) == NULL){  /* K.ycomplex */
    nycomplex = 0;
  }
  else{
    nycomplex = mxGetM(K_FIELD) * mxGetN(K_FIELD);
    ycomplexPr = mxGetPr(K_FIELD);
  }
/* --------------------------------------------------
   GET STATISTICS:
   -------------------------------------------------- */
  nblk = 1 + cK.lorN + cK.sdpN;
/* ------------------------------------------------------------
   Allocate working array blkstart(nblk+1).
   ------------------------------------------------------------ */
  blkstart = (int *) mxCalloc(nblk + 1, sizeof(int));
  getblkstart(blkstart, cK);
/* ------------------------------------------------------------
   Create output structure K
   ------------------------------------------------------------ */
  K_OUT = mxCreateStructMatrix(1, 1, NKFIELDS, Kfieldnames);
/* ------------------------------------------------------------
   Create fields K.{l,q,s, rsdpN, blkstart}
   ------------------------------------------------------------ */
  setScalarField(K_OUT,"l",cK.lpN);                          /* K.l */
  K_FIELD = mxCreateDoubleMatrix(1,cK.lorN, mxREAL);         /* K.q */
  memcpy(mxGetPr(K_FIELD), cK.lorNL, cK.lorN * sizeof(double));
  mxSetField(K_OUT,0,"q",K_FIELD);
  K_FIELD = mxCreateDoubleMatrix(1,cK.sdpN, mxREAL);         /* K.s */
  memcpy(mxGetPr(K_FIELD), cK.sdpNL, cK.sdpN * sizeof(double));
  mxSetField(K_OUT,0,"s",K_FIELD);
  setScalarField(K_OUT,"rsdpN",cK.rsdpN);                    /* K.rsdpN */
  K_FIELD = mxCreateDoubleMatrix(nblk+1,1, mxREAL);          /* K.blkstart */
  blkstartPr = mxGetPr(K_FIELD);
  for(i = 0; i <= nblk; i++)
    blkstartPr[i] = blkstart[i] + 1;
  mxSetField(K_OUT,0,"blkstart",K_FIELD);
  K_FIELD = mxCreateDoubleMatrix(1,nycomplex, mxREAL);     /* K.ycomplex */
  memcpy(mxGetPr(K_FIELD), ycomplexPr, nycomplex * sizeof(double));
  mxSetField(K_OUT,0,"ycomplex",K_FIELD);
/* ------------------------------------------------------------
   Create fields K.{rLen, hLen, qMaxn, rMaxN, hMaxn}
   ------------------------------------------------------------ */
  setScalarField(K_OUT,"rLen",cK.rLen);                          /* K.rLen */
  setScalarField(K_OUT,"hLen",cK.hLen);                          /* K.hLen */
  setScalarField(K_OUT,"qMaxn",cK.qMaxn);                        /* K.qMaxn */
  setScalarField(K_OUT,"rMaxn",cK.rMaxn);                        /* K.rMaxn */
  setScalarField(K_OUT,"hMaxn",cK.hMaxn);                        /* K.hMaxn */
/* ------------------------------------------------------------
   Release working array
   ------------------------------------------------------------ */
  mxFree(blkstart);
}
