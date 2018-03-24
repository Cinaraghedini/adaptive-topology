/*
   [vfrm.q,dq,ud.qdet] = qframev(x,ux.qdet,dIN,udIN.qdet,dscl,vfrm.lab,...
     sqrt(detv),K);
   Computes 1. dq = D(d) dscl with d in K, Lorentz part.
            2. vfrm = v2/(sqrt(2)*norm(v2)), v is spec geom mean, Lorentz.

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
#include <math.h>
#include "mex.h"
#include "triuaux.h"
#include "blksdp.h"

#define VFRM_OUT myplhs[0]
#define DQ_OUT myplhs[1]
#define QDET_OUT myplhs[2]
#define NPAROUT 3

#define X_IN prhs[0]
#define RDETX_IN prhs[1]
#define D_IN prhs[2]
#define RDETD_IN prhs[3]
#define DSCL_IN prhs[4]
#define VLAB_IN prhs[5]
#define RDETV_IN prhs[6]
#define K_IN prhs[7]
#define NPARIN 8

/* ************************************************************
   PROCEDURE myqlmul : LORENTZ SCALE z = D(x)y (full version)
     z=D(x)y = [x'*y / sqrt(2);  alpha * x(2:n) + rdetx * y(2:n)],
     where alpha = (z1+rdetx*y1) / (x(1)+ sqrt(2) * rdetx)
   INPUT
     x,y - full n x 1
     rdetx - sqrt(det(x))
     n - order of x,y,z.
     denom - x[0] + M_SQRT2 * rdetx.
   OUTPUT
     z - full n x 1. Let z := D(x)y.
   ************************************************************ */
void myqlmul(double *z,const double *x,const double *y,
             const double rdetx,const int n, const double denom)
{
  double z1;
/* ------------------------------------------------------------
   z1 = x'*y / sqrt(2),
   alpha = (z1+rdetx*y1) / denom
   z(1) = z1, z(2:n) = alpha * x(2:n) + rdetx * y(2:n).
    ------------------------------------------------------------ */
  z1 = (z[0] = realdot(x,y,n) / M_SQRT2);
  scalarmul(z+1,(z1 + rdetx * y[0]) / denom,x+1,n-1);
  addscalarmul(z+1,rdetx,y+1,n-1);
}

/* ************************************************************
   TIME-CRITICAL PROCEDURE -- isscalardiv
   Computes  x /= alpha using loop-unrolling.
   INPUT
     alpha - scalar
     n - order
   UPDATED
     x - Length n, on output, x /= alpha.
   ************************************************************ */
void isscalardiv(double *x, const double alpha,const int n)
{
 int k;

 for(k = 0; k < n-3; k++){                 /* LEVEL 4 */
   x[k] /= alpha; k++;
   x[k] /= alpha; k++;
   x[k] /= alpha; k++;
   x[k] /= alpha;
 }
 if(k < n-1){                              /* LEVEL 2 */
   x[k] /= alpha; k++;
   x[k] /= alpha; k++;
 }
 if(k < n)                                 /* LEVEL 1 */
   x[k] /= alpha;
}

/* ************************************************************
   PROCEDURE qframev - Computes the frame of v = D(dnew)\(D(d)x)
     avoiding cancelation from ill-conditioning of d.
   INPUT
     xscl  - primal iterate in scaled space, length n.
     dscl, rdetdscl - new scaling in scaled space.
     d     - length n, original scaling.
     dnew1, rdetdnew - dnew = D(d)*dscl. dnew1 = dnew[0].
     denom - d[0] + sqrt(2) * rdetd
     vlab - length 2, eigK(v), vlab(1) <= vlab(2).
     n  - order of Lorentz block.
   OUTPUT
     vfrm - length n-1 vector. On output, vfrm = v(2:n)/(vlab2-vlab1)
   ************************************************************ */
void qframev(double *vfrm, const double *xscl, const double *dscl,
             const double rdetdscl, const double *d,
             const double dnew1, const double rdetdnew,
             const double denom, const double *vlab, const int n)
{
  double denomnew, x1,v1, del1,del2,del3, alpha1,alpha2;

  if(vlab[0] >= vlab[1])
    fzeros(vfrm, n-1);
  else{
    denomnew = dnew1 + M_SQRT2 * rdetdnew;
    x1 = realdot(d,xscl,n) / M_SQRT2;              /* x = D(d)*xscl */
    v1 = (vlab[0] + vlab[1]) / M_SQRT2;
/* ------------------------------------------------------------
   v is a linear combination of xscl, d and dscl. Compute the multipliers.
   ------------------------------------------------------------ */
    del1 = M_SQRT2 * rdetdscl - dscl[0];        /* 0 if d is identity */
    del2 = xscl[0] - rdetdscl * v1;                /* 0 if dscl is identity */
    del3 = M_SQRT2 * xscl[0] - dscl[0] * v1;       /* idem */
    alpha1 = (x1 + rdetdnew * v1) / denomnew;
    alpha2 = ((x1 * del1) + (dnew1 * del2 + rdetdnew * del3)) /
      (denom * denomnew);
/* ------------------------------------------------------------
   rdetdscl * v2 = x2 - alpha1 * dscl2 + alpha2 * d2.
   vfrm = v2 / (vlab2 - vlab1).
   ------------------------------------------------------------ */
    memcpy(vfrm, xscl+1, (n-1) * sizeof(double));
    subscalarmul(vfrm, alpha1, dscl+1, n-1);
    addscalarmul(vfrm, alpha2, d+1,n-1);
    isscalardiv(vfrm, rdetdscl * (vlab[1]-vlab[0]), n-1);
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
  mxArray *myplhs[NPAROUT];
  int i,k, nk;
  double *vfrm, *dq, *qdet, *fwork;
  const double *d,*rdetd,*x,*rdetx, *dscl, *vlab, *rdetv;
  double denom;
  coneK cK;
/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  if(nrhs < NPARIN)
    mexErrMsgTxt("qframev requires more input arguments.");
  if (nlhs > NPAROUT)
    mexErrMsgTxt("qframev generates less output arguments.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Get scale data: (d,rdetd) and inputs x, dscl, vlab.
   ------------------------------------------------------------ */
  if(mxGetM(X_IN) * mxGetN(X_IN) < cK.lpN + cK.qDim)            /* x */
    mexErrMsgTxt("x size mismatch");
  x = mxGetPr(X_IN) + cK.lpN;              /* skip LP part */
  if(mxGetM(RDETX_IN) * mxGetN(RDETX_IN) != cK.lorN)            /* rdetx */
    mexErrMsgTxt("rdetx size mismatch");
  rdetx = mxGetPr(RDETX_IN);
  if(mxGetM(D_IN) * mxGetN(D_IN) < cK.lpN + cK.qDim)            /* d */
    mexErrMsgTxt("d size mismatch");
  d = mxGetPr(D_IN) + cK.lpN;              /* skip LP part */
  if(mxGetM(RDETD_IN) * mxGetN(RDETD_IN) != cK.lorN)            /* rdetd */
    mexErrMsgTxt("rdetd size mismatch");
  rdetd = mxGetPr(RDETD_IN);
  if(mxGetM(DSCL_IN) * mxGetN(DSCL_IN) != cK.qDim)              /* dscl */
    mexErrMsgTxt("dscl size mismatch");
  dscl = mxGetPr(DSCL_IN);
  if(mxGetM(VLAB_IN) * mxGetN(VLAB_IN) < cK.lpN + cK.lorN)      /* vlab */
    mexErrMsgTxt("vlab size mismatch");
  vlab = mxGetPr(VLAB_IN) + cK.lpN;              /* skip LP part */
  if(mxGetM(RDETV_IN) * mxGetN(RDETV_IN) != cK.lorN)      /* rdetv */
    mexErrMsgTxt("rdetv size mismatch");
  rdetv = mxGetPr(RDETV_IN);
/* ------------------------------------------------------------
   Allocate outputs VFRM_OUT(qDim-lorN), DQ_OUT(qDim), qdet(lorN)
   ------------------------------------------------------------ */
  VFRM_OUT =  mxCreateDoubleMatrix(cK.qDim - cK.lorN, 1, mxREAL);
  vfrm = mxGetPr(VFRM_OUT);
  DQ_OUT =  mxCreateDoubleMatrix(cK.qDim, 1, mxREAL);
  dq = mxGetPr(DQ_OUT);
  QDET_OUT = mxCreateDoubleMatrix(cK.lorN, 1, mxREAL);
  qdet = mxGetPr(QDET_OUT);
/* ------------------------------------------------------------
   Allocate working arrays: double fwork(lorN).
   ------------------------------------------------------------ */
  fwork = (double *) mxCalloc(MAX(1,cK.lorN), sizeof(double));
/* ------------------------------------------------------------
   The actual job is done here:
   fwork = "rdetdscl" = rdet(x) ./ rdetv,
   qdet = rdetd .* rdetdscl,
   ------------------------------------------------------------ */
  realHadadiv(fwork, rdetx, rdetv, cK.lorN);
  realHadamard(qdet, rdetd, fwork, cK.lorN);
  for(k = 0; k < cK.lorN; k++){               /* LORENTZ */
    nk = cK.lorNL[k];
/* ------------------------------------------------------------
   Let dq = D(d)*dscl
   ------------------------------------------------------------ */
    denom = d[0] + M_SQRT2 * rdetd[k];
    myqlmul(dq, d,dscl,rdetd[k],nk,denom);
/* ------------------------------------------------------------
   Let vfrm = frame(v)
   ------------------------------------------------------------ */
    qframev(vfrm, x, dscl, fwork[k], d, dq[0], qdet[k], denom, vlab, nk);
    d += nk; dscl += nk; dq += nk;
    x += nk; vlab += 2; vfrm += nk-1;
  }
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
  mxFree(fwork);
/* ------------------------------------------------------------
   Copy requested output parameters (at least 1), release others.
   ------------------------------------------------------------ */
  i = MAX(nlhs, 1);
  memcpy(plhs,myplhs, i * sizeof(mxArray *));
  for(; i < NPAROUT; i++)
    mxDestroyArray(myplhs[i]);
}
