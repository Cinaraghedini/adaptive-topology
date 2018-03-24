/*
    z = invjmulK(xlab,xfrm, y, K)
    solves x jmul z = y, with x = XFRM*xlab in spectral factors.

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

#include <string.h>
#include "mex.h"
#include "blksdp.h"
#include "reflect.h"

/* z = invjmulK(xlab,xfrm, y, K) */
#define Z_OUT plhs[0]
#define NPAROUT 1

#define X_IN prhs[0]
#define FRM_IN prhs[1]
#define Y_IN prhs[2]
#define K_IN prhs[3]
#define NPARIN 4

/* ============================================================
   LORENTZ: SOLVE {Z : X = Y JMUL Z}
   ============================================================ */

/* ************************************************************
   PROCEDURE qjdiv : LORENTZ INVERSE JMUL
   INPUT
     x - full n x 1
     y1, y2 - spectral values of y
     n - order of x,z.
     f - length n-1 frame of y.
   OUTPUT
     z - full n x 1. Solves y jmul z = x.
   ************************************************************ */
void qjdiv(double *z,const double *x,const double y1,const double y2,
           const int n,const double *frame)
{
 double fx2, r2y1, difyr2, z1;
/* ------------------------------------------------------------
   fx2 = f'*x(2:n);
   r2y1 = 2 / (y(1) +y(2));
   difyr2 = (y(2) - y(1)) / sqrt(2);
   ------------------------------------------------------------ */
 fx2 = realdot(frame,x+1,n-1);
 r2y1 = 2 / (y1 +y2);
 difyr2 = (y2 - y1) / M_SQRT2;
/* ------------------------------------------------------------
   z = [(x(1)/r2y1 - difyr2*fx2) / (y1*y2); ...
   r2y1 * x(2:n) - r2y1 * z(1)*difyr2*f]
   ------------------------------------------------------------ */
 z1 = (x[0] / r2y1 - difyr2 * fx2) / (y1 * y2);
 z[0] = z1;
 scalarmul(z+1, r2y1,x+1,n-1);
 addscalarmul(z+1,- r2y1 * z1 * difyr2,frame,n-1);
}

/* ============================================================ 
   PSD: z(i.j) = 2*y(i,j)/(xi + xj)   (Z and Y are matrices, x is vector)
   ============================================================ */

/* ************************************************************
   PROCEDURE diagjdiv  -   yij := 2*y(i,j)/(xi + xj) for i >= j.
     The strict upper triangular of Y is left unchanged.
   INPUT
     x - length n vector
     n - order of x and square y matrix
   UPDATED
     y - full n x n, on return tril(yNEW) = 2*y(i,j)/(xi + xj), i >= j.
   ************************************************************ */
void diagjdiv(double *y,const double *x,const int n)
{
  int i,j;
  double xj;

/* ------------------------------------------------------------
   For j=0..n-1:
   for i=j..n-1:  let y(i,j) *= 2/(xi + xj)
   ------------------------------------------------------------ */
 for(j = 0; j < n; y+=n, j++){
   xj = x[j];
   y[j] /= xj;
   for(i = j+1; i < n; i++)
     y[i] *= 2 / (x[i] + xj);
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
 const mxArray *FRM_FIELD;
 int k, nk, nksqr, lenfull, lendiag, sdpdim, qsize;
 double *z, *fwork;
 const double *x,*y, *frmq, *frms, *beta;
 coneK cK;

/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
 if(nrhs < NPARIN)
   mexErrMsgTxt("invjmulK requires more input arguments.");
 if(nlhs > NPAROUT)
   mexErrMsgTxt("invjmulK generates 1 output argument.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
 conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Get statistics of cone K structure
   ------------------------------------------------------------ */
 sdpdim = cK.rDim + cK.hDim;
 qsize = sdpdim + cK.hLen;
 lenfull = cK.lpN +  cK.qDim + sdpdim;
 lendiag = cK.lpN + 2 * cK.lorN + cK.rLen + cK.hLen;
/* ------------------------------------------------------------
   Get inputs x, frm, y.
   ------------------------------------------------------------ */
 if(mxGetM(Y_IN) * mxGetN(Y_IN) != lenfull)        /* validate x, y */
   mexErrMsgTxt("size y mismatch.");
 if(mxGetM(X_IN) * mxGetN(X_IN) != lendiag)
   mexErrMsgTxt("size xlab mismatch.");
 if(mxIsSparse(X_IN) | mxIsSparse(Y_IN))
   mexErrMsgTxt("Sparse inputs not supported by this version of invjmulK.");
 x = mxGetPr(X_IN);                                /* get x and y */
 y = mxGetPr(Y_IN);
 if(!mxIsStruct(FRM_IN))
   mexErrMsgTxt("Parameter `xfrm' should be a structure.");
 if( (FRM_FIELD = mxGetField(FRM_IN,0,"q")) == NULL)         /* xfrm.q */
   mexErrMsgTxt("Field xfrm.q missing.");
 frmq = mxGetPr(FRM_FIELD);
 if( (FRM_FIELD = mxGetField(FRM_IN,0,"s")) == NULL)         /* xfrm.s */
   mexErrMsgTxt("Field xfrm.s missing.");
 if(mxGetM(FRM_FIELD) * mxGetN(FRM_FIELD) != qsize)
   mexErrMsgTxt("size xfrm.s mismatch.");
 frms = mxGetPr(FRM_FIELD);
/* ------------------------------------------------------------
   Allocate output Z
   ------------------------------------------------------------ */
 Z_OUT =  mxCreateDoubleMatrix(lenfull, 1, mxREAL);
 z = mxGetPr(Z_OUT);
/* ------------------------------------------------------------
   Allocate working array fwork(max(rmaxn,2*hmaxn))
   ------------------------------------------------------------ */
  fwork = (double *) mxCalloc(MAX(1,MAX(cK.rMaxn,2*cK.hMaxn)),sizeof(double));
/* ------------------------------------------------------------
   The actual job is done here:.
   ------------------------------------------------------------ */
/* ------------------------------------------------------------
   LP: z = y./x
   ------------------------------------------------------------ */
 realHadadiv(z, y,x,cK.lpN);
 z += cK.lpN;             /* Next, point to lorentz & sdp blocks */
 x += cK.lpN; y += cK.lpN;
/* ------------------------------------------------------------
   LORENTZ: solve  z s.t. x jmul z = y
   ------------------------------------------------------------ */
 for(k = 0; k < cK.lorN; k++){
   nk = cK.lorNL[k];
   qjdiv(z,y,x[0],x[1],nk, frmq);
   z += nk; y += nk; x += 2;
   frmq += (nk-1);
 }
/* ------------------------------------------------------------
   PSD: Since X = Q'*diag(x)*Q, we have XZ+ZX = 2Y iff
    qzqt(i,j) = 2*qyqt(i,j)/(xi + xj),  qzqt = Q*z*Q', qyqt = Q*y*Q'.
   ------------------------------------------------------------ */
 for(k = 0; k < cK.rsdpN; k++){
   nk = cK.sdpNL[k];
   nksqr = SQR(nk);
/* ------------------------------------------------------------
   Let z = Q*y*Q'
   ------------------------------------------------------------ */
   memcpy(z, y, nksqr * sizeof(double));
   beta = frms + nksqr - nk;
   qxqt(z, beta, frms, nk, fwork);
/* ------------------------------------------------------------
   Solve diag(x) jmul zNEW = zOLD
   ------------------------------------------------------------ */
   diagjdiv(z,x,nk);
/* ------------------------------------------------------------
   Let zFINAL = Q'*z*Q  (back into old eig-basis)
   ------------------------------------------------------------ */
   qtxq(z, beta, frms, nk, fwork);
   tril2sym(z,nk);
   z += nksqr; y += nksqr; frms += nksqr;
   x += nk;
 }
 for(; k < cK.sdpN; k++){                    /* complex Hermitian */
   nk = cK.sdpNL[k];
   nksqr = SQR(nk);
/* ------------------------------------------------------------
   Let z = Q*y*Q'
   ------------------------------------------------------------ */
   memcpy(z, y, 2 * nksqr * sizeof(double));
   beta = frms + 2 * nksqr;                   /* beta = frms(:,2*n+1) */
   prpiqxqt(z,z+nksqr, beta, frms,frms+nksqr, nk, fwork);
/* ------------------------------------------------------------
   Solve diag(x) jmul zNEW = zOLD. Since x is real, we can handle
   the real and imaginary parts seperately.
   ------------------------------------------------------------ */
   diagjdiv(z,x,nk);            /* real part */
   diagjdiv(z+nksqr,x,nk);      /* imaginary part */
/* ------------------------------------------------------------
   Let zFINAL = Q'*z*Q  (back into old eig-basis)
   ------------------------------------------------------------ */
   prpiqtxq(z,z+nksqr, beta, frms,frms+nksqr, nk, fwork);
   tril2herm(z,z+nksqr,nk);
   nksqr += nksqr;
   z += nksqr; y += nksqr; frms += nksqr + nk;          /* skip also beta. */
   x += nk;
 }
/* ------------------------------------------------------------
   Release working array
   ------------------------------------------------------------ */
 mxFree(fwork);
}
