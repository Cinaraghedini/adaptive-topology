#include "givens.h"

/* ************************************************************
   PROCEDURE givensrot - apply sequence of givens rotations to
   a vector.
   INPUT
     g - length n: each entry is a givens rotation [x, y;y,-x], x^2+y^2=1.
       We apply first g[n-1] to z[n-1:n], up to g[0] to z[0:1].
     n - order of g, i.e. number of givens rotations.
   UPDATED
     z - length n+1 vector, to be rotated n times.
   ************************************************************ */
void givensrot(double *z, const twodouble *g, const int n)
{
  twodouble gi;
  double z1, z2;
  int i;

  z2 = z[n];
  for(i = n-1; i >= 0; i--){
    gi = g[i];
    z1 = z[i];
/* ------------------------------------------------------------
   [z1NEW;    [x,  y;       [z1;
    z2NEW] :=  y, -x]   *    z2]
   ------------------------------------------------------------ */
    z[i+1] = gi.y * z1 - gi.x * z2;           /* z2NEW */
    z2 = gi.x * z1 + gi.y * z2;               /* z1NEW, is z2 in iter --i */
  }
  z[0] = z2;
}

/* ************************************************************
   PROCEDURE prpigivensrot - apply sequence of givens rotations to
   a vector. Complex case.
   INPUT
     g - length n: each entry is a givens rotation [conj(x), y;y,-x],
       |x|^2+y^2=1. y is real.
       We apply first g[n-1] to z[n-1:n], up to g[0] to z[0:1].
     n - order of g, i.e. number of givens rotations.
   UPDATED
     z,zpi - length n+1 vector, to be rotated n times.
   ************************************************************ */
void prpigivensrot(double *z,double *zpi, const tridouble *g, const int n)
{
  tridouble gi;
  double z1, z2, z1im,z2im;
  int i;

  z2 = z[n];
  z2im = zpi[n];
  for(i = n-1; i >= 0; i--){
    gi = g[i];
    z1 = z[i];
    z1im = zpi[i];
/* ------------------------------------------------------------
   [z1NEW;    [conj(x),  y;       [z1;
    z2NEW] :=        y, -x]   *    z2]
   ------------------------------------------------------------ */
/* z2NEW = y*z1 - x*z2 , y is real. */
    z[i+1] = gi.y * z1 - gi.x * z2 + gi.xim * z2im;
    zpi[i+1] = gi.y * z1im - gi.x * z2im - gi.xim * z2;
/* z1NEW, is z2 in iter --i. z1NEW = conj(x)*z1 + y*z2, y is real. */
    z2 = gi.x * z1 + gi.xim * z1im + gi.y * z2;
    z2im = gi.x * z1im - gi.xim * z1 + gi.y * z2im;
  }
  z[0] = z2;
  zpi[0] = z2im;
}

/* ************************************************************
   PROCEDURE givensrotuj - apply sequence of givens rotations to
   a vector, whose last affected entry is now 0. Typical for
   re-inserting columns in a U-factor. Same as "givensrot", except that
   z[n]=0 by assumption on input.
   INPUT
     g - length n: each entry is a givens rotation [x, y;y,-x], x^2+y^2=1.
       We apply first g[n-1] to z[n-1:n], up to g[0] to z[0:1].
     n - order of g, i.e. number of givens rotations.
   UPDATED
     z - length n+1 vector, to be rotated n times. On input, z[n] = 0
      by assumption (actual contents irrelevant).
   ************************************************************ */
void givensrotuj(double *z, const twodouble *g, const int n)
{
  twodouble gi;
  double z1, z2;
  int i;

  if(n < 1)
    return;
/* ------------------------------------------------------------
   [z2;      [x,  y;       [z[n-1];
    z[n]] :=  y, -x]   *    0]      = z[n-1] * [x;y]
   ------------------------------------------------------------ */
  z2 = z[n-1];
  z[n] = z2 * ((g+n-1)->y);
  z2 *= (g+n-1)->x;
  for(i = n-2; i >= 0; i--){
    gi = g[i];
    z1 = z[i];
/* ------------------------------------------------------------
   [z1NEW;    [x,  y;       [z1;
    z2NEW] :=  y, -x]   *    z2]
   ------------------------------------------------------------ */
    z[i+1] = gi.y * z1 - gi.x * z2;           /* z2NEW */
    z2 = gi.x * z1 + gi.y * z2;               /* z1NEW, is z2 in iter --i */
  }
  z[0] = z2;
}

/* ************************************************************
   PROCEDURE prpigivensrotuj - apply sequence of givens rotations to
   a vector, whose last affected entry is now 0, and the preceding
   entry is real. On output, the last effected entry will be real.
   Typical for re-inserting columns in a U-factor. Same as "givensrot",
   except that z[n-1] is real and z[n]=0 by assumption on input.
   INPUT
     g - length n: each entry is a givens rotation [conj(x), y;y,-x],
       |x|^2+y^2=1, y is real.
       We apply first g[n-1] to z[n-1:n], up to g[0] to z[0:1].
     n - order of g, i.e. number of givens rotations.
   UPDATED
     z - length n+1 vector, to be rotated n times. On input,
      {zpi[n-1],z[n],zpi[n]} = 0 by assumption (actual contents irrelevant).
   ************************************************************ */
void prpigivensrotuj(double *z,double *zpi, const tridouble *g, const int n)
{
  tridouble gi;
  double z1, z2, z1im,z2im;
  int i;

  if(n < 1)
    return;
/* ------------------------------------------------------------
   [z2;      [conj(x),  y;       [z[n-1];
    z[n]] :=        y, -x]   *    0]      = z[n-1] * [conj(x);y],
   where z[n-1] is real.
   ------------------------------------------------------------ */
  z2 = z[n-1];
  z[n] = z2 * ((g+n-1)->y);
  z2im = -z2 * (g+n-1)->xim;                /* z[n-1] * conj(x) */
  z2 *= (g+n-1)->x;
  for(i = n-2; i >= 0; i--){
    gi = g[i];
    z1 = z[i];
    z1im = zpi[i];
/* ------------------------------------------------------------
   [z1NEW;    [conj(x),  y;       [z1;
    z2NEW] :=        y, -x]   *    z2]
   ------------------------------------------------------------ */
/* z2NEW = y*z1 - x*z2 , y is real. */
    z[i+1] = gi.y * z1 - gi.x * z2 + gi.xim * z2im;
    zpi[i+1] = gi.y * z1im - gi.x * z2im - gi.xim * z2;
/* z1NEW, is z2 in iter --i. z1NEW = conj(x)*z1 + y*z2, y is real. */
    z2 = gi.x * z1 + gi.xim * z1im + gi.y * z2;
    z2im = gi.x * z1im - gi.xim * z1 + gi.y * z2im;
  }
  z[0] = z2;
  zpi[0] = z2im;
}
