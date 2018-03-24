/* ************************************************************
%                                        Ablkjc = partitA(At,blkstart)
% PARTITA  Partition columns of A according to the subscripts listed in blkstart.
%
% SEE ALSO sedumi
% ******************** INTERNAL FUNCTION OF SEDUMI ********************
function Ablkjc = partitA(At,blkstart) --  Partition columns of A
  according to the subscripts listed in blkstart


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
#include <math.h>      /* floor and log */
#include "mex.h"
#include "blksdp.h"

#define ABLKJC_OUT plhs[0]
#define NPAROUT 1

#define AT_IN prhs[0]
#define BLKSTART_IN prhs[1]
#define NPARIN 2

void intadd(int *x, const int y, const int n)
{
  int i;
  for(i = 0; i < n; i++)
    x[i] += y;
}

/* ************************************************************
   PROCEDURE partitA
   INPUT
     Ajc, Air, m  - sparse N x m matrix
     blkstart, nblk - length nblk integer array of subscripts.
     iwsize - length of iwork, iwsize = floor(log(1+nblk)/log(2)).
   OUTPUT
     Ablkjc - length (nblk+2)*m array. Rows 1+(1:nblk) list 1st nonzero
       with subscript at or beyond blkstart.
   WORK
     cfound - length nblk char work array
     iwork  - length iwsize = floor(log(1+nblk)/log(2)) work array.
   ************************************************************ */
void partitA(int *Ablkjc, const int *Ajc,const int *Air,
             const int *blkstart, const int m,const int nblk,
             const int iwsize, char *cfound, int *iwork)
{
  int j, L;
  L = nblk+2;
  for(j = 0; j < m; j++)
    intmbsearch(Ablkjc + j*L, cfound, Air+Ajc[j], Ajc[j+1]-Ajc[j],
		blkstart, nblk, iwork, iwsize);
  for(j = 0; j < m; j++)
    intadd(Ablkjc + j*L, Ajc[j],L);
}

/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
                 const int nrhs, const mxArray *prhs[])
{
  jcir At;
  int i,j, nblk,m, L, iwsize;
  int *iwork, *Ablkjc, *blkstart;
  const int *rowj;
  double *AblkjcPr;
  const double *blkstartPr;
  char *cwork;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  if(nrhs < NPARIN)
    mexErrMsgTxt("partitA requires more input arguments.");
  if (nlhs > NPAROUT)
    mexErrMsgTxt("partitA produces less output arguments.");
/* --------------------------------------------------
   GET inputs At, blkstart
   -------------------------------------------------- */
  if(!mxIsSparse(AT_IN))
    mexErrMsgTxt("At must be a sparse matrix.");
  At.jc = mxGetJc(AT_IN);
  At.ir = mxGetIr(AT_IN);
  m = mxGetN(AT_IN);
  nblk = mxGetM(BLKSTART_IN) * mxGetN(BLKSTART_IN);
  blkstartPr = mxGetPr(BLKSTART_IN);
/* ------------------------------------------------------------
   Allocate working array Ablkjc((nblk+2) * m), iwork(log_2(1+nblk)),
   blkstart(nblk)
   ------------------------------------------------------------ */
  iwsize = floor(log(1+nblk)/log(2));
  iwork = (int *) mxCalloc(MAX(iwsize,1), sizeof(int));
  Ablkjc = (int *) mxCalloc(MAX((nblk+2)*m,1), sizeof(int));
  blkstart = (int *) mxCalloc(MAX(nblk,1), sizeof(int));
  cwork = (char *) mxCalloc(MAX(nblk,1), sizeof(char));
/* ------------------------------------------------------------
   Translate blkstart from Fortran-double to C-int
   ------------------------------------------------------------ */
  for(i = 0; i < nblk; i++){                         /* to integers */
    j = blkstartPr[i];
    blkstart[i] = --j;
  }
/* ------------------------------------------------------------
   The real job:
   ------------------------------------------------------------ */
  partitA(Ablkjc, At.jc,At.ir, blkstart, m,nblk, iwsize,cwork,iwork);
/* ------------------------------------------------------------
   Create output Ablkjc m x nblk.
   ------------------------------------------------------------ */
  ABLKJC_OUT = mxCreateDoubleMatrix(m, nblk, mxREAL);
  AblkjcPr = mxGetPr(ABLKJC_OUT);
  rowj = Ablkjc;
  L = nblk+2;
  for(j = 0; j < nblk; j++){
    ++rowj;
    for(i = 0; i < m; i++)
      AblkjcPr[i] = rowj[i*L];      /* convert int to double */
    AblkjcPr += m;
  }
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
  mxFree(cwork);
  mxFree(iwork);
  mxFree(Ablkjc);
  mxFree(blkstart);
}
