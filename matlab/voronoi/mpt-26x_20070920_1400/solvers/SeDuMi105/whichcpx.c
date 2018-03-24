/*
%                                                    cpx = whichcpx(K)
% WHICHCPX  yields structure cpx.{f,q,r,x}
%
% SEE ALSO sedumi
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

#include "mex.h"
#include "blksdp.h"

#define CPX_OUT plhs[0]
#define K_IN prhs[0]


/* ************************************************************
   PROCEDURE whichcpx - Split xcomplex into [ifr, xcomplex], and
     determine number of extra entries for each q- and r-cone.
   INPUT
     nxcomplex - length(xcomplex on input)
     frlpN, lorN, rconeN - K.f+K.l, length(K.q), length(K.r), resp.
   UPDATED
     xcomplex - Length nxcomplex integer array. On input, the indices
       tbat are to be interpreted as complex entries. On output, it has
       only those that are part of norm-bounded q/r-cone entries.
       Length nxcomplex - n.
     lorNL - Length lorN array, listing order of each q-cone on input.
       On output, it has the extra entries per q-cone from the IM-parts.
     rconeNL - Length lorN array, listing order of each r-cone on input.
       On output, it has the extra entries per r-cone from the IM-parts.
   OUTPUT
     ifr - length n <= nxcomplex array. Lists the x_i's of which
       the IM(x_i)-part should be treated as free variable.
   RETURNS n=number of free imaginary parts, length(ifr).
   ************************************************************ */
int whichcpx(int *ifr, int *xcomplex, int *lorNL, int *rconeNL,
	     const int nxcomplex, const int frlpN, const int lorN,
	     const int rconeN)
{
  int n, ix, i,j,k,lastj, ixOld;
  if(nxcomplex <= 0)
    return 0;
  n = 0; ix = 0;   /* target index into ifr, xcomplex */
/* ------------------------------------------------------------
   Complex free and LP-nonneg variables result in free IM-parts
   ------------------------------------------------------------ */
  for(i = 0; i < nxcomplex; i++){        /* i is source into xcomplex */
    if( (j = xcomplex[i]) >= frlpN)
      break;
    ifr[n++] = j;
  }
/* ------------------------------------------------------------
   For the Lorentz q-vectors, the 1st entry results in free IM-part,
   but the remaining in extra q-entries for that cone.
   ------------------------------------------------------------ */
  lastj = frlpN;
  for(k = 0; k < lorN; k++){
    if(j == lastj){
      ifr[n++] = j;                 /* complex 1st entry */
      i++;
    }
    lastj += lorNL[k];              /* point beyond current block */
    ixOld = ix;
    for(; i < nxcomplex; i++){
      if( (j = xcomplex[i]) >= lastj)
	break;
      xcomplex[ix++] = j;           /* complex within kth Lorentz block */
    }
    lorNL[k] = ix - ixOld;           /* number of extra entries in q-cone */
  }
/* ------------------------------------------------------------
   For the Lorentz r-vectors, the 1st 2 entries
   can be complex; it there results in a free IM-part.
   Complex values for the remaining entries results in extra r-entries
   for that cone.
   ------------------------------------------------------------ */
  for(k = 0; k < rconeN; k++){
    if(j == lastj || j == lastj + 1){
      ifr[n++] = j;                 /* complex 1st entry */
      if(++i < nxcomplex)
	if(xcomplex[i] == lastj + 1){
	  ifr[n++] = j;
	  i++;
	}
    }
    lastj += rconeNL[k];
    ixOld = ix;
    for(; i < nxcomplex; i++){
      if( (j = xcomplex[i]) >= lastj)
	break;
      xcomplex[ix++] = j;
    }
    rconeNL[k] = ix - ixOld;        /* number of extra entries in r-cone */
  }
/* ------------------------------------------------------------
   Success: return number of extra free variables (from IM-parts)
   ------------------------------------------------------------ */
  return n;
}


/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
#define NCPX_FIELDS 4
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   x = whichcpx(K)
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
  int i,j,iwsiz, nxcomplex, cpxf;
  int *iwork, *lorNL, *rconeNL, *xcomplex;
  double *myPr;
  const double *xcomplexPr;
  mxArray *MY_FIELD;
  const char *CPXFieldnames[] = {"f", "q", "r", "x"};
  coneK cK;
/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  if(nrhs < 1)
    mexErrMsgTxt("whichcpx requires 1 input argument.");
  if(nlhs > 1)
    mexErrMsgTxt("whichcpx generates 1 output argument.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
  if( (MY_FIELD = mxGetField(K_IN,0,"xcomplex")) == NULL){  /* K.xcomplex */
    nxcomplex = 0;
  }
  else{
    nxcomplex = mxGetM(MY_FIELD) * mxGetN(MY_FIELD);
    xcomplexPr = mxGetPr(MY_FIELD);
  }
  if(nxcomplex > 0){
/* ------------------------------------------------------------
   ALLOCATE working arrays:
   iwork(2*nxcomplex+lorN+rconeN))
   ------------------------------------------------------------ */
    iwsiz = 2* nxcomplex + cK.lorN + cK.rconeN;
    iwork = (int *) mxCalloc(MAX(iwsiz,1), sizeof(int));
    xcomplex = iwork + nxcomplex;
    lorNL = xcomplex + nxcomplex;
    rconeNL = lorNL + cK.lorN;
/* ------------------------------------------------------------
   Convert double to int
   ------------------------------------------------------------ */
    for(i = 0; i < nxcomplex; i++){
      j = xcomplexPr[i];                       /* double to int */
      xcomplex[i] = --j;                       /* Fortran to C */
    }
    for(i = 0; i < cK.lorN; i++)
      lorNL[i] = cK.lorNL[i];                  /* double to int */
    for(i = 0; i < cK.rconeN; i++)
      rconeNL[i] = cK.rconeNL[i];              /* double to int */
/* ------------------------------------------------------------
   The real work:
   ------------------------------------------------------------ */
    cpxf = whichcpx(iwork, xcomplex, lorNL,rconeNL,
		    nxcomplex, cK.frN + cK.lpN, cK.lorN, cK.rconeN);
    nxcomplex -= cpxf;
  }
/* ------------------------------------------------------------
   If xcomplex = []:
   ------------------------------------------------------------ */
  else{
    cpxf = 0;
    iwsiz = 0;
  }
/* ------------------------------------------------------------
   Create output structure CPX
   ------------------------------------------------------------ */
  CPX_OUT = mxCreateStructMatrix(1, 1, NCPX_FIELDS, CPXFieldnames);
  MY_FIELD = mxCreateDoubleMatrix(cpxf,1,mxREAL);      /* cpx.f */
  myPr = mxGetPr(MY_FIELD);
  for(i = 0; i < cpxf; i++)
    myPr[i] = 1.0 + iwork[i];                          /* int to double */
  mxSetField(CPX_OUT, 0,"f", MY_FIELD);
  MY_FIELD = mxCreateDoubleMatrix(cK.lorN,1,mxREAL);      /* cpx.q */
  if(iwsiz > 0){
    myPr = mxGetPr(MY_FIELD);
    for(i = 0; i < cK.lorN; i++)
      myPr[i] = lorNL[i];                              /* int to double */
  }
  mxSetField(CPX_OUT, 0,"q", MY_FIELD);
  MY_FIELD = mxCreateDoubleMatrix(cK.rconeN,1,mxREAL);      /* cpx.r */
  if(iwsiz > 0){
    myPr = mxGetPr(MY_FIELD);
    for(i = 0; i < cK.rconeN; i++)
      myPr[i] = rconeNL[i];                              /* int to double */
  }
  mxSetField(CPX_OUT, 0,"r", MY_FIELD);
  MY_FIELD = mxCreateDoubleMatrix(nxcomplex,1,mxREAL);      /* cpx.x */
  myPr = mxGetPr(MY_FIELD);
  for(i = 0; i < nxcomplex; i++)
    myPr[i] = 1.0 + xcomplex[i];                       /* int to double */
  mxSetField(CPX_OUT, 0,"x", MY_FIELD);
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
  if(iwsiz > 0)
    mxFree(iwork);
}
