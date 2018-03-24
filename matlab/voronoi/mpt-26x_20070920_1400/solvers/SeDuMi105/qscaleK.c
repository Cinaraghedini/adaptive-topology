/*
   y = qscaleK(d,rdetd,x,K)
   Computes y = D(d) x with d in K, Lorentz part.

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
#include "triuaux.h"
#include "blksdp.h"

/*    y = scaleK(d,ud,x,K) */

#define Y_OUT plhs[0]
#define NPAROUT 1

#define D_IN prhs[0]
#define RDETD_IN prhs[1]
#define X_IN prhs[2]
#define K_IN prhs[3]
#define NPARIN 4

/* ============================================================
   LORENTZ OPERATIONS
   ============================================================ */
/* ************************************************************
   PROCEDURE qscale : LORENTZ SCALE D(x)y = z + mu * x (full version)
     mu = (y1+alpha)/sqrt(2), z = rdetx * [alpha; y(2:n)],
     where alpha = (x(2:n)'*y(2:n)) / (x(1)+ sqrt(2) * rdetx)
   INPUT
     x,y - full n x 1
     rdetx - sqrt(det(x))
     n - order of x,y,z.
   OUTPUT
     z - full n x 1. Let z := rdetx * [alpha; y(2:n)].
   RETURNS
     mu = (y1+alpha)/sqrt(2).
   ************************************************************ */
double qscale(double *z,const double *x,const double *y,
              const double rdetx,const int n)
{
 double alpha, mu;
/* ------------------------------------------------------------
   alpha = (x(2:n)'*y(2:n)) / (x(1)+ sqrt(2) * rdetx)
   ------------------------------------------------------------ */
 alpha = realdot(x+1,y+1,n-1) / (x[0] + M_SQRT2 * rdetx);
/* ------------------------------------------------------------
   z = rdetx * [alpha; y(2:n)].
   ------------------------------------------------------------ */
 z[0] = rdetx * alpha;
 scalarmul(z+1,rdetx,y+1,n-1);
/* ------------------------------------------------------------
   RETURN mu = (y1+alpha)/sqrt(2).
   ------------------------------------------------------------ */
 return (y[0] + alpha) / M_SQRT2;
}

/* ************************************************************
   PROCEDURE qlmul : LORENTZ SCALE z = D(x)y (full version)
     z=D(x)y = [x'*y / sqrt(2);  mu * x(2:n) + rdetx * y(2:n)],
     where mu = (z(1)+rdetx*y1) / (x(1)+ sqrt(2) * rdetx)
   INPUT
     x,y - full n x 1
     rdetx - sqrt(det(x))
     n - order of x,y,z.
   OUTPUT
     z - full n x 1. Let z := D(x)y.
   ************************************************************ */
void qlmul(double *z,const double *x,const double *y,
	   const double rdetx,const int n)
{
  double mu;
  mu = qscale(z, x,y,rdetx,n);
  addscalarmul(z,mu,x,n);
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
  int i,k, nk;
  double *y;
  const double *d,*rdetd,*x;
  coneK cK;
/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  if(nrhs < NPARIN)
    mexErrMsgTxt("qscaleK requires more input arguments.");
  if (nlhs > NPAROUT)
    mexErrMsgTxt("qscaleK generates 1 output argument.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Get scale data: (d,rdetd) and input x.
   ------------------------------------------------------------ */
  if(mxGetM(D_IN) * mxGetN(D_IN) < cK.lpN + cK.qDim)
    mexErrMsgTxt("d size mismatch");
  d = mxGetPr(D_IN) + cK.lpN;              /* skip LP part */
  if(mxGetM(RDETD_IN) * mxGetN(RDETD_IN) != cK.lorN)
    mexErrMsgTxt("rdetx size mismatch");
  rdetd = mxGetPr(RDETD_IN);
  if(mxGetM(X_IN) * mxGetN(X_IN) != cK.qDim)
    mexErrMsgTxt("x size mismatch");
  x = mxGetPr(X_IN);
/* ------------------------------------------------------------
   Allocate output Y
   ------------------------------------------------------------ */
  Y_OUT =  mxCreateDoubleMatrix(cK.qDim, 1, mxREAL);
  y = mxGetPr(Y_OUT);
/* ------------------------------------------------------------
   The actual job is done here: y=D(d)x, Lorentz part.
   ------------------------------------------------------------ */
  for(k = 0; k < cK.lorN; k++){               /* LORENTZ */
    nk = cK.lorNL[k];
    qlmul(y, d,x,rdetd[k],nk);
    y += nk; x += nk; d += nk;
  }
}
