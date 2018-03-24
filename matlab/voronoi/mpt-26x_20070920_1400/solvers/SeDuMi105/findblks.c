/* ************************************************************
%                         Ablk = findblks(At,Ablkjc,blk0,blk1,blkstart)
% FINDBLKS  Find nonzero blocks
%   in A, with subscripts per column bounded bij Ablkjc([blk0,blk1]),
%   block partitioned by blkstart.
%   If blk0 < 1 (blk1 > size(Ablkjc,2)) then start (stop) searching at column
%   start (end) of A.
%
%  SEE ALSO partitA.
% ******************** INTERNAL FUNCTION OF SEDUMI ********************

function Ablk = findblks(At,Ablkjc,blk0,blk1,blkstart) --  Find nonzero blocks

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
#include <math.h>
#include "mex.h"
#include "blksdp.h"

#define ABLK_OUT plhs[0]
#define NPAROUT 1

#define AT_IN prhs[0]
#define ABLKJC_IN prhs[1]
#define BLK0_IN prhs[2]
#define BLK1_IN prhs[3]
#define BLKSTART_IN prhs[4]
#define NPARIN 5

/* ************************************************************
   PROCEDURE findblks
   INPUT
     Ajc1, Ajc2, Air, m  - sparse N x m matrix, we only consider nonzeros
       in Air[Ajc1[k]:Ajc2[k]-1], k=0:m-1.
     blkstart, nblk - length nblk integer array of subscripts.
     blkstartm1 - length nblk array, blkstartm1[k]=blkstart[k]-1
     iwsize - length of iwork, iwsize = nblk+2+floor(log(1+nblk)/log(2)).
   OUTPUT
     Ablkjc, Ablkir - sparse nblk x m matrix, less than sum(Ajc2-Ajc1)
       nonzeros.
   WORK
     cfound - length nblk char work array
     iwork  - length iwsize = nblk+2+floor(log(1+nblk)/log(2)) work array.
   ************************************************************ */
void findblks(int *Ablkir, int *Ablkjc, const int *Ajc1,const int *Ajc2,
              const int *Air, const int *blkstart, const int *blkstartm1,
              const int m,const int nblk,
              int iwsize, char *cfound, int *iwork)
{
  int i,j,inz,ajnnz;
  int *ipos;
/* ------------------------------------------------------------
   Partition working array into ipos(nblk+2), iwork.
   ------------------------------------------------------------ */
  ipos = iwork;
  iwork += nblk+2;
  iwsize -= nblk+2;
  inz = 0;
  for(j = 0; j < m; j++){
    Ablkjc[j] = inz;
/* ------------------------------------------------------------
   If A(:,j) has more nonzeros than blkstart, we search blkstart
   ------------------------------------------------------------ */
    if((ajnnz = Ajc2[j]-Ajc1[j]) > nblk){
      intmbsearch(ipos, cfound, Air+Ajc1[j], ajnnz,
                  blkstart, nblk, iwork, iwsize);
      for(i = 0; i < nblk; i++)
        if(ipos[i+2] > ipos[i+1])
          Ablkir[inz++] = i;
    }
    else{
/* ------------------------------------------------------------
   If A(:,j) has less nonzeros than blkstart, we search those nonzeros
   within blkstartm1. The position of the nonzero is then the block number+1.
   ------------------------------------------------------------ */
      intmbsearch(ipos, cfound, blkstartm1,nblk, Air+Ajc1[j],ajnnz,
                  iwork, iwsize);
      for(i = 0; i < ajnnz; i++)
        if(ipos[i+1] > ipos[i])         /* New block number ? */
          Ablkir[inz++] = ipos[i+1] - 1;   /* ipos is block number + 1 */
    }
  }
/* ------------------------------------------------------------
   Close last column of Ablk
   ------------------------------------------------------------ */
  Ablkjc[m] = inz;
}

/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
                 const int nrhs, const mxArray *prhs[])
{
  jcir At, Ablk;
  int i,j, nblk,m, blknnz, njc, iwsize, blk0,blk1;
  int *iwork, *Ajc, *blkstart;
  const double *blkstartPr, *AjcPr;
  char *cwork;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  if(nrhs < NPARIN)
    mexErrMsgTxt("findblks requires more input arguments.");
  if (nlhs > NPAROUT)
    mexErrMsgTxt("findblks produces less output arguments.");
/* --------------------------------------------------
   GET inputs At, blkstart, Ablkjc, blk0, blk1
   -------------------------------------------------- */
  if(!mxIsSparse(AT_IN))
    mexErrMsgTxt("At must be a sparse matrix.");
  At.jc = mxGetJc(AT_IN);
  At.ir = mxGetIr(AT_IN);
  m = mxGetN(AT_IN);
  nblk = mxGetM(BLKSTART_IN) * mxGetN(BLKSTART_IN) - 1;
  blkstartPr = mxGetPr(BLKSTART_IN);
  AjcPr = mxGetPr(ABLKJC_IN);
  if(m != mxGetM(ABLKJC_IN))
    mexErrMsgTxt("Ablkjc size mismatch.");
  njc = mxGetN(ABLKJC_IN);
  blk0 = mxGetScalar(BLK0_IN);           /* double to int */
  --blk0;                                /* Fortran to C */
  if(mxGetM(BLK1_IN) * mxGetN(BLK1_IN) != 1)
    blk1 = njc;                           /*default to end */
  else{
    blk1 = mxGetScalar(BLK1_IN);   /* double to int (thus inf not allowed) */
    --blk1;                                /* Fortran to C */
  }
/* ------------------------------------------------------------
   Allocate working array iwork(nblk+2+log_2(1+nblk)),
   blkstart(2*nblk), Ajc(2*m)
   char cwork(nblk)
   ------------------------------------------------------------ */
  iwsize = nblk + 2 + floor(log(1+nblk)/log(2));
  iwork = (int *) mxCalloc(iwsize, sizeof(int));
  blkstart = (int *) mxCalloc(MAX(2*nblk,1), sizeof(int));
  Ajc = (int *) mxCalloc(MAX(2*m,1), sizeof(int));
  cwork = (char *) mxCalloc(MAX(nblk,1), sizeof(char));
/* ------------------------------------------------------------
   Translate blkstart from Fortran-double to C-int
   ------------------------------------------------------------ */
  for(i = 0; i < nblk; i++){                         /* to integers */
    j = blkstartPr[i];
    blkstart[i] = --j;
    blkstart[nblk+i] = --j;          /* blkstart minus 1 */
  }
/* ------------------------------------------------------------
   Convert Ajc from double to int:
   ------------------------------------------------------------ */
  if(blk0 < 0)
    memcpy(Ajc,At.jc,m*sizeof(int));          /* default: start of column */
  else if(blk0 >= njc)
    mexErrMsgTxt("Ablkjc size mismatches blk0.");
  else
    for(i = 0; i < m; i++){                         /* to integers */
      Ajc[i] = AjcPr[m*blk0 + i];
    }
  if(blk1 >= njc)
    memcpy(Ajc+m,At.jc+1,m*sizeof(int));      /* default: end of column */
  else if(blk1 < 0)
    mexErrMsgTxt("blk1 must be positive.");
  else
    for(i = 0; i < m; i++){                         /* to integers */
      Ajc[m+i] = AjcPr[blk1*m + i];
    }
/* ------------------------------------------------------------
   Ablk = sparse(nblk,m,blknnz);
   ------------------------------------------------------------ */
  blknnz = 0;
  for(i = 0; i < m; i++)
    blknnz += Ajc[m+i]-Ajc[i];              /* upper bound on nnz blocks */
  blknnz = MAX(blknnz,1);
  ABLK_OUT = mxCreateSparse(nblk,m, blknnz,mxREAL);
  Ablk.jc = mxGetJc(ABLK_OUT);
  Ablk.ir = mxGetIr(ABLK_OUT);
/* ------------------------------------------------------------
   The real job:
   ------------------------------------------------------------ */
  findblks(Ablk.ir,Ablk.jc, Ajc,Ajc+m,At.ir, blkstart,blkstart+nblk,
           m,nblk, iwsize,cwork,iwork);
/* ------------------------------------------------------------
   REALLOC (shrink) Ablk to Ablk.jc[m] nonzeros.
   ------------------------------------------------------------ */
  mxAssert(Ablk.jc[m] <= blknnz,"");
  blknnz = MAX(Ablk.jc[m],1);
  if((Ablk.ir = (int *) mxRealloc(Ablk.ir, blknnz * sizeof(int))) == NULL)
    mexErrMsgTxt("Memory allocation error");
  if((Ablk.pr = (double *) mxRealloc(mxGetPr(ABLK_OUT), blknnz*sizeof(double)))
     == NULL)
    mexErrMsgTxt("Memory allocation error");
  mxSetPr(ABLK_OUT,Ablk.pr);
  mxSetIr(ABLK_OUT,Ablk.ir);
  mxSetNzmax(ABLK_OUT,blknnz);
  for(i = 0; i < blknnz; i++)
    Ablk.pr[i] = 1.0;
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
  mxFree(cwork);
  mxFree(iwork);
  mxFree(Ajc);
  mxFree(blkstart);
}
