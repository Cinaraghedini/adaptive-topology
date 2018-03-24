/*
   y = jmulK(x,y, K)  for full x,y, or
   z = jmulK(x,y, K,f) with x in spectral space and Lorentz frame f.
   z = x jmul y (Jordan multiplication).
   Argument "f" only needed if x is in spectral space.

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

#define Z_OUT plhs[0]
#define X_IN prhs[0]
#define Y_IN prhs[1]
#define K_IN prhs[2]
#define F_IN prhs[3]


/* ============================================================
   LORENTZ: Z = X JMUL Y
   ============================================================ */

/* ************************************************************
   PROCEDURE qjmul : LORENTZ JMUL
   INPUT
     x, y - full n x 1
     n - order of x,y,z.
   OUTPUT
     z - full n x 1. Let z := x jmul y.
   ************************************************************ */
void qjmul(double *z,const double *x,const double *y, const int n)
{
  /* ------------------------------------------------------------
     z = [x'*y; (x(1)*y(2:n) + y(1)*x(2:n))] / sqrt(2);
     ------------------------------------------------------------ */
  z[0] = realdot(x,y,n) / M_SQRT2;
  ++z;
  scalarmul(z, x[0] / M_SQRT2,y+1,n-1);
  addscalarmul(z, y[0] / M_SQRT2,x+1,n-1);
}

/* ************************************************************
   PROCEDURE frqjmul : LORENTZ SCALE (framed version)
   INPUT
     x - full n x 1
     y1, y2 - spectral values of y
     n - order of x,z.
     frame - length n-1 frame of y
   OUTPUT
     z - full n x 1. Let z := x jmul y.
   ************************************************************ */
void frqjmul(double *z,const double *x,const double y1,const double y2,
	     const int n,const double *frame)
{
 double fx2, y1r2, difyr2, x1;
  /* ------------------------------------------------------------
     fx2 = f'*x(2:n);
     y1r2 = (y(1) +y(2)) / 2;
     difyr2 = (y(2) - y(1)) / sqrt(2);
     ------------------------------------------------------------ */
 x1 = x[0];
 x++;
 fx2 = realdot(frame,x,n-1);
 y1r2 = (y1 +y2) / 2;
 difyr2 = (y2 - y1) / M_SQRT2;
 /* ------------------------------------------------------------
    z = [(y1r2 * x(1) + difyr2*fx2);  (y1r2 * x(2:n) + (x(1)*difyr2)*f)]
    ------------------------------------------------------------ */
 z[0] = y1r2 * x1 + difyr2 * fx2;
 ++z;
 scalarmul(z,y1r2,x,n-1);
 addscalarmul(z,x1 * difyr2,frame,n-1);
}

/* ============================================================ 
   SDP: Z = tril(X*Y+Y*X)/2
   ============================================================ */

/* ************************************************************
   PROCEDURE symjmul(z,x,y,n)  --
     Z = tril(X * Y + Y * X) / 2, with X,Y symmetric.
     The strict upper triangular of Z is undefined.
   INPUT
     x,y - symmetric n x n
     n - order of square x,y,z matrices
   UPDATED
     z - full n x n, on entry arbitrary, on return tril(Z) = tril(XY+YX)/2.
   ************************************************************ */
void symjmul(double *z,const double *x,const double *y,const int n)
{
 int i,j,icol;

 /* ------------------------------------------------------------
    For i=0..n-1:
        for j=i..n-1:  let z(j,i) = (x(:,i)'*y(:,j) + y(:,i)'*x(:,j))/2
    ------------------------------------------------------------ */
 for(j = 0; j < n; z += n, x+=n, y+=n, j++){
   z[j] = realdot(x,y,n);
   for(i = j+1, icol = n; i < n; icol+=n, i++)
     z[i] = (realdot(x+icol,y,n) + realdot(x,y+icol,n)) / 2;
 }
}

/* ************************************************************
   PROCEDURES hermjmul --
     Z = tril(X * Y + Y * X) / 2, with X,Y Hermitian. X = [RE X,IM X].
     Only stict lower triangular and real diagonal of R are defined.
   INPUT
     x,y - Hermitian; full 2*(n x n).
     n - order of square x,y,z matrices
   UPDATED
     z - full 2*(n x n), on entry arbitrary, on return
       tril(Z) = tril(XY + YX)/2.
   ************************************************************ */
void hermjmul(double *z,const double *x,const double *y,const int n)
{
 int i,j,icol,jcol,nsqr;
 const double *xpi,*ypi;

 nsqr = SQR(n);
 xpi = x + nsqr; ypi = y + nsqr;
 /* ------------------------------------------------------------
    For i=0..n-1:
        for j=i..n-1:  let z(j,i) = (x(:,i)'*y(:,j) + y(:,i)'*x(:,j))/2
    ------------------------------------------------------------ */
 for(j = 0, jcol=0; j < n; z += n, jcol+=n, j++){
   z[j] = realdot(x+jcol,y+jcol,n) + realdot(xpi+jcol,ypi+jcol,n);
   for(i = j+1, icol = jcol; i < n; i++){
     icol += n;
     z[i] = (realdot(x+icol,y+jcol,n) + realdot(x+jcol,y+icol,n)
       + realdot(xpi+icol,ypi+jcol,n) + realdot(xpi+jcol,ypi+icol,n)) / 2;
   }
 }
 for(j = 0, jcol=0; j < n; z += n, jcol+=n, j++){
   for(i = j+1, icol = jcol; i < n; i++){
     icol += n;
     z[i] = (realdot(x+icol,ypi+jcol,n)
	     - realdot(xpi+icol,y+jcol,n)
	     - realdot(x+jcol,ypi+icol,n)
	     + realdot(xpi+jcol,y+icol,n)) / 2;
   }
 }
}

/* ************************************************************
   PROCEDURE diagjmul(z,x,y,n)  --
     Z = tril(X * Y + Y * X) / 2, with Y symmetric and X diagonal, i.e.
       zij = (xi + xj) * yij /2
     The strict upper triangular of Z is undefined.
   INPUT
     x - length n vector
     y - symmetric n x n
     n - order of square y,z matrices
   UPDATED
     z - full n x n, on entry arbitrary, on return tril(Z) = tril(XY+YX)/2.
   ************************************************************ */
void diagjmul(double *z,const double *x,const double *y,const int n)
{
 int i,j;
 double xj;

 /* ------------------------------------------------------------
    For j=0..n-1:
        for i=j..n-1:  let z(i,j) = (xi + xj) * y(i,j) / 2
    ------------------------------------------------------------ */
 for(j = 0; j < n; z += n, y+=n, j++){
   xj = x[j];
   z[j] = xj * y[j];
   for(i = j+1; i < n; i++)
     z[i] = (x[i] + xj) * y[i] / 2;
 }
}

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
    z = jmulK(x,y, K)  for full x,y, or
    z = jmulK(x,y, K,f) with x in spectral space and Lorentz frame f.
    z = x jmul y (Jordan multiplication).
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
 int k, nk, nksqr, lenfull, lendiag, lenx;
 double *z;
 const double *x,*y, *frame;
 coneK cK;

/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
 if(nrhs < 3)
   mexErrMsgTxt("jmulK requires at least 3 input arguments.");
 if(nlhs > 1)
   mexErrMsgTxt("jmulK generates 1 output argument.");
 /* ------------------------------------------------------------
    Disassemble cone K structure
    ------------------------------------------------------------ */
 conepars(K_IN, &cK);
 /* ------------------------------------------------------------
    Get statistics of cone K structure
    ------------------------------------------------------------ */
 lenfull = cK.lpN +  cK.qDim + cK.rDim + cK.hDim;
 lendiag = cK.lpN + 2 * cK.lorN + cK.rLen + cK.hLen;
 /* ------------------------------------------------------------
    Get inputs x and y.
    ------------------------------------------------------------ */
 x = mxGetPr(X_IN);
 y = mxGetPr(Y_IN);
 if(mxIsSparse(X_IN) | mxIsSparse(Y_IN))
   mexErrMsgTxt("Sparse inputs not supported by this version of jmulK.");
 /* ------------------------------------------------------------
    Validate input:
    length(y) == lenfull ?
    length(x) in {lenfull, lendiag} ? (lendiag -> get frame)
    NB: the idea is that X,Y are Hermitian, but this we don't check.
    ------------------------------------------------------------ */
 if(mxGetM(Y_IN) * mxGetN(Y_IN) != lenfull)
   mexErrMsgTxt("size y mismatch.");
 lenx = mxGetM(X_IN) * mxGetN(X_IN);
 if(lenx != lenfull){
   if(lenx != lendiag)
     mexErrMsgTxt("size x mismatch.");
   else if(cK.lorN){
     if(nrhs < 4)
        mexErrMsgTxt("Input argument f required.");
     else{
       frame = mxGetPr(F_IN);
     }
   }
 }
 /* ------------------------------------------------------------
    Allocate output Z
    ------------------------------------------------------------ */
 Z_OUT =  mxCreateDoubleMatrix(lenfull, 1, mxREAL);
 z = mxGetPr(Z_OUT);
 /* ------------------------------------------------------------
    The actual job is done here:.
    ------------------------------------------------------------ */
 if(cK.lpN){
   /* ------------------------------------------------------------
      LP: z = x.*y
      ------------------------------------------------------------ */
   realHadamard(z, x,y,cK.lpN);
   z += cK.lpN;             /* Next, point to lorentz & sdp blocks */
   x += cK.lpN; y += cK.lpN;
 }
 if(cK.lorN){
   /* ------------------------------------------------------------
      LORENTZ: z = x jmul y
      (I) x in R^N full:
      ------------------------------------------------------------ */
   if(lenx > lendiag)
     for(k = 0; k < cK.lorN; k++){
       nk = cK.lorNL[k];
       qjmul(z,x,y,nk);
       z += nk; x += nk; y += nk;
     }
   else
     /* ------------------------------------------------------------
	(II) x in R^n spectral: using Lorentz frame
	------------------------------------------------------------ */
     for(k = 0; k < cK.lorN; k++){
       nk = cK.lorNL[k];
       frqjmul(z,y,x[0],x[1],nk, frame);
       z += nk; y += nk; x += 2;
       frame += (nk-1);
     }
 }
 /* ------------------------------------------------------------
    PSD: Z = (XY + YX)/2
    (I) full
    ------------------------------------------------------------ */
 if(lenx > lendiag){
   for(k=0; k < cK.rsdpN; k++){                /* real symmetric */
     nk = cK.sdpNL[k];
     symjmul(z,x,y,nk);
     tril2sym(z,nk);
     nksqr = SQR(nk);
     z += nksqr; x += nksqr;
     y += nksqr;
   }
   for(; k < cK.sdpN; k++){                    /* complex Hermitian */
     nk = cK.sdpNL[k];
     nksqr = SQR(nk);
     hermjmul(z,x,y,nk);
     tril2herm(z,z+nksqr,nk);
     nksqr += nksqr;
     z += nksqr; x += nksqr;
     y += nksqr;
   }
 }
 else{
   /* ------------------------------------------------------------
      (II) z(i.j) = y(i,j)*(xi + xj)/2 if x spectral.
      ------------------------------------------------------------ */
   for(k=0; k < cK.rsdpN; k++){
     nk = cK.sdpNL[k];
     diagjmul(z,x,y,nk);
     tril2sym(z,nk);
     nksqr = SQR(nk);
     z += nksqr; y += nksqr;
     x += nk;
   }
   for(; k < cK.sdpN; k++){
     nk = cK.sdpNL[k];
     diagjmul(z,x,y,nk);
     nksqr = SQR(nk);
     diagjmul(z+nksqr,x,y+nksqr,nk);
     tril2herm(z,z+nksqr,nk);
     nksqr += nksqr;
     z += nksqr; y += nksqr;
     x += nk;
   }
 }
}
