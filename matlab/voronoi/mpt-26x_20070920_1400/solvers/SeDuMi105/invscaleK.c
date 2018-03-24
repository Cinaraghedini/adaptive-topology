/*
   y = invscaleK(d,ud,x,K,transp)
   Computes y = D(d^{-1}) x with d in K.
   !transp then Y = Ld' \ X / Ld
   transp == 1 then Y = Ud' \ X / Ud

   Uses pivot ordering ud.perm if available and nonempty.

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

/*    y = invscaleK(d,ud,x,K,transp) */

#define Y_OUT plhs[0]
#define NPAROUT 1

#define D_IN prhs[0]
#define UD_IN prhs[1]
#define X_IN prhs[2]
#define K_IN prhs[3]
#define NPARINMIN 4
#define TRANSP_IN prhs[4]
#define NPARIN 5

/* ************************************************************
   PROCEDURE qldiv : LORENTZ SCALE z = D(x)\y (full version)
    D(x)\y = (1/det x) * [x'Jy/sqrt(2); rdetx * y2-alpha*x2],
    where alpha = (x'Jy/sqrt(2) + rdetx*y1) / (x(1)+ sqrt(2) * rdetx)
   INPUT
     x,y - full n x 1
     rdetx - sqrt(det(x))
     n - order of x,y,z.
   OUTPUT
     z - full n x 1. Let z := D(x)^{-1}y.
   ************************************************************ */
void qldiv(double *z,const double *x,const double *y,
	   const double rdetx,const int n)
{
 double alpha,x1,y1,z1;
/* ------------------------------------------------------------
   z1 = x'*J*y / (sqrt(2) * det x),
   alpha = (z1+y1/rdetx) / (x(1)+ sqrt(2) * rdetx)
   ------------------------------------------------------------ */
 x1 = x[0]; y1 = y[0];
 z1 = (x1*y1 - realdot(x+1,y+1,n-1)) / (M_SQRT2 * SQR(rdetx));
 alpha = (z1 + y1 / rdetx) / (x1 + M_SQRT2 * rdetx);
/* ------------------------------------------------------------
   z(1) = z1, z(2:n) = y(2:n)/rdetx - alpha * x(2:n).
   ------------------------------------------------------------ */
 z[0] = z1;
 scalarmul(z+1,-alpha,x+1,n-1);
 addscalarmul(z+1,1/rdetx,y+1,n-1);
}

/* ************************************************************
   PROCEDURE invscaleK - Computes y = D(d^{-1})x.
   For PSD, uses D=U'*U factorization, (transp == 0) Y = U\X/U'
   or (transp == 1) Y = U'\X/U.
   INPUT
     x - length N(K) input vector.
     d - scaling vector, only LP and Lorentz part needed.
     ud - Cholesky factor of d for PSD part (after PERM ordering).
     qdetd - sqrt(det(d)) for Lorentz part.
     perm - ordering: UD=chol(d(perm,perm)), for numerical stability.
     cK   - structure describing symmetric cone K.
   OUTPUT
     y - length N(K) output vector, y=D(d)x.
   WORK
     fwork - length 2 * max(rmaxn^2,2*hmaxn^2) working vector.
   ************************************************************ */
void invscaleK(double *y, const double *d, const double *ud,
               const double *qdetd, const int *perm, const double *x,
               const coneK cK, const char transp, double *fwork)
{
  int k,nk,nksqr;
  double *z,*zpi;
  char use_pivot;
/* ------------------------------------------------------------
   Partition fwork into fwork(psdblk) and z(psdblk), where
   psdblk = max(rmaxn^2,2*hmaxn^2).
   Let zpi = z + hmaxn^2.
   ------------------------------------------------------------ */
  use_pivot = (perm != (const int *) NULL);
  z = fwork + MAX(SQR(cK.rMaxn),2*SQR(cK.hMaxn));
  zpi = z + SQR(cK.hMaxn);
/* ------------------------------------------------------------
   LP: y = x ./ d
   ------------------------------------------------------------ */
  realHadadiv(y, x,d,cK.lpN);
  y += cK.lpN;             /* Next, point to lorentz & sdp blocks */
  d += cK.lpN; x += cK.lpN;
/* ------------------------------------------------------------
   LORENTZ: y = D(d^{-1})* x
   ------------------------------------------------------------ */
  for(k = 0; k < cK.lorN; k++){
    nk = cK.lorNL[k];
    qldiv(y, d,x,qdetd[k],nk);
    y += nk; d += nk; x +=nk;
  }
/* ------------------------------------------------------------
   PSD, !transp: triu(Y) = triu(Ld' \ (X(perm,perm) / Ld)).
   Needs ony tril(X/Ld).
   ------------------------------------------------------------ */
  if(!transp){
    if(use_pivot){            /* with pivoting */
      for(k = 0; k < cK.rsdpN; k++){                /* real symmetric */
        nk = cK.sdpNL[k];
        invltxl(z,ud,x,nk,fwork);
        tril2sym(z,nk);
        invmatperm(y,z,perm,nk);            /* Y(perm,perm) = Z */
        nksqr = SQR(nk);
        y += nksqr; ud += nksqr;
        x += nksqr; perm += nk;
      }
      for(; k < cK.sdpN; k++){                    /* complex Hermitian */
        nk = cK.sdpNL[k];
        nksqr = SQR(nk);
        prpiinvltxl(z,zpi,ud,ud+nksqr,x,x+nksqr,nk,fwork);
        tril2herm(z,zpi,nk);
        invmatperm(y,z,perm,nk);            /* Y(perm,perm) = Z */
        invmatperm(y+nksqr,zpi,perm,nk);    /* imaginary part */
        nksqr += nksqr;
        y += nksqr; ud += nksqr;
        x += nksqr; perm += nk;
      }
    }
    else{ /* without pivoting */
      for(k = 0; k < cK.rsdpN; k++){                /* real symmetric */
        nk = cK.sdpNL[k];
        invltxl(y,ud,x,nk,fwork);
        tril2sym(y,nk);
        nksqr = SQR(nk);
        y += nksqr; ud += nksqr;
        x += nksqr;
      }
      for(; k < cK.sdpN; k++){                    /* complex Hermitian */
        nk = cK.sdpNL[k];
        nksqr = SQR(nk);
        prpiinvltxl(y,y+nksqr,ud,ud+nksqr,x,x+nksqr,nk,fwork);
        tril2herm(y,y+nksqr,nk);
        nksqr += nksqr;
        y += nksqr; ud += nksqr;
        x += nksqr;
      }
    } /* no pivot */
  } /* !transp */
  else{
/* ------------------------------------------------------------
   PSD, transp: triu(Y) = triu(Ud' \ (X(perm,perm) / Ud)).
   Needs ony triu(X/Ud).
   ------------------------------------------------------------ */
    if(use_pivot){            /* with pivoting */
      for(k = 0; k < cK.rsdpN; k++){                /* real symmetric */
        nk = cK.sdpNL[k];
        matperm(z,x,perm,nk);                    /* z = x(perm,perm) */
        invutxu(y,ud,z,nk,fwork);
        triu2sym(y,nk);
        nksqr = SQR(nk);
        y += nksqr; ud += nksqr;
        x += nksqr; perm += nk;
      }
      for(; k < cK.sdpN; k++){                    /* complex Hermitian */
        nk = cK.sdpNL[k];
        nksqr = SQR(nk);
        matperm(z,x,perm,nk);                    /* z = x(perm,perm) */
        matperm(zpi,x+nksqr,perm,nk);                    /* imaginary part */
        prpiinvutxu(y,y+nksqr,ud,ud+nksqr,z,zpi,nk,fwork);
        triu2herm(y,y+nksqr,nk);
        nksqr += nksqr;
        y += nksqr; ud += nksqr;
        x += nksqr; perm += nk;
      }
    }
    else{ /* without pivoting */
      for(k = 0; k < cK.rsdpN; k++){                /* real symmetric */
        nk = cK.sdpNL[k];
        invutxu(y,ud,x,nk,fwork);
        triu2sym(y,nk);
        nksqr = SQR(nk);
        y += nksqr; ud += nksqr;
        x += nksqr;
      }
      for(; k < cK.sdpN; k++){                    /* complex Hermitian */
        nk = cK.sdpNL[k];
        nksqr = SQR(nk);
        prpiinvutxu(y,y+nksqr,ud,ud+nksqr,x,x+nksqr,nk,fwork);
        triu2herm(y,y+nksqr,nk);
        nksqr += nksqr;
        y += nksqr; ud += nksqr;
        x += nksqr;
      }
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
  const mxArray *UD_FIELD;
  int lenfull, lend, lenud, sdplen, transp, fwsiz, i,k;
  double *fwork, *y,*permPr;
  const double *d,*x,*ud, *qdetd;
  int *perm, *iwork;
  coneK cK;
  char use_pivot;

/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  if(nrhs < NPARIN){
    transp = 0;
    if(nrhs < NPARINMIN)
      mexErrMsgTxt("invscaleK requires more input arguments.");
  }
  else
    transp = mxGetScalar(TRANSP_IN);
  if (nlhs > NPAROUT)
    mexErrMsgTxt("invscaleK generates 1 output argument.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Get statistics of cone K structure
   ------------------------------------------------------------ */
  lenfull = cK.lpN +  cK.qDim + cK.rDim + cK.hDim;
  lend = cK.lpN +  cK.qDim;
  lenud = cK.rDim + cK.hDim;
  sdplen = cK.rLen + cK.hLen;
/* ------------------------------------------------------------
   Get scale data: (d,ud) and input x.
   ------------------------------------------------------------ */
  d = mxGetPr(D_IN);
  if(!mxIsStruct(UD_IN))
    mexErrMsgTxt("Parameter `ud' should be a structure.");
  if( (UD_FIELD = mxGetField(UD_IN,0,"u")) == NULL)         /* ud.u */
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
  if( (UD_FIELD = mxGetField(UD_IN,0,"qdet")) == NULL)         /* ud.qdet */
    mexErrMsgTxt("Field ud.qdet missing.");
  if(mxGetM(UD_FIELD) * mxGetN(UD_FIELD) != cK.lorN)
    mexErrMsgTxt("ud.qdet size mismatch.");
  qdetd = mxGetPr(UD_FIELD);
/* ------------------------------------------------------------
   Get input x
   ------------------------------------------------------------ */
  if(mxIsSparse(X_IN))
    mexErrMsgTxt("Sparse x not supported by this version of scaleK.");
  x = mxGetPr(X_IN);
/* ------------------------------------------------------------
   Validate input:
   length(d) == lend ?  length(x) == lenfull ?
   NB: the idea is that X is Hermitian, but this we don't check.
   ------------------------------------------------------------ */
  if(mxGetM(D_IN) * mxGetN(D_IN) < lend)
    mexErrMsgTxt("size d mismatch.");
  if(mxGetM(X_IN) * mxGetN(X_IN) != lenfull)
    mexErrMsgTxt("size x mismatch.");
/* ------------------------------------------------------------
   Allocate output Y
   ------------------------------------------------------------ */
  Y_OUT =  mxCreateDoubleMatrix(lenfull, 1, mxREAL);
  y = mxGetPr(Y_OUT);
/* ------------------------------------------------------------
   Allocate fwork 2 * [ max(cK.rMaxn^2, 2*cK.hMaxn^2) ]
   iwork = int(rLen+hLen)
   ------------------------------------------------------------ */
  fwsiz = MAX(SQR(cK.rMaxn),2*SQR(cK.hMaxn));
  fwork = (double *) mxCalloc( MAX(1,2 * fwsiz), sizeof(double));
  iwork = (int *) mxCalloc( MAX(1, cK.rLen + cK.hLen), sizeof(int) );
/* ------------------------------------------------------------
   Convert Fortran to C-style in perm:
   ------------------------------------------------------------ */
  if(use_pivot){
    perm = iwork;
    for(k = 0; k < cK.rLen + cK.hLen; k++){
      i = permPr[k];
      perm[k] = --i;
    }
  }
  else
    perm = (int *) NULL;
/* ------------------------------------------------------------
   The actual job is done here:.
   ------------------------------------------------------------ */
  invscaleK(y, d,ud,qdetd,perm,x,cK,transp, fwork);
/* ------------------------------------------------------------
   Release working arrays.
   ------------------------------------------------------------ */
  mxFree(fwork);
  mxFree(iwork);
}
