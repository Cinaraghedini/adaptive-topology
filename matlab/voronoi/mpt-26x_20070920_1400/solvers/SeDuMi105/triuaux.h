#include "blksdp.h"
#ifndef TRIUAUX
#define TRIUAUX

void matperm(double *y, const double *x, const int *perm, const int n);
void invmatperm(double *y, const double *x, const int *perm, const int n);
void realltxl(double *y, const double *l, const double *x,
              const int n, double *xl);
void prpiltxl(double *y,double *ypi, const double *ld,const double *ldpi,
              const double *x,const double *xpi,
              const int n, double *xld);
void realutxu(double *y, const double *u, const double *x,
              const int n, double *xu);
void prpiutxu(double *y, double *ypi, const double *u, const double *upi,
              const double *x, const double *xpi,
              const int n, double *xu);
void invutxu(double *y, const double *u, const double *x,
              const int m, double *xu);
void prpiinvutxu(double *y,double *ypi, const double *u,const double *upi,
             const double *x,const double *xpi, const int m, double *xu);
void invltxl(double *y, const double *l, const double *x,
              const int m, double *xl);
void prpiinvltxl(double *y,double *ypi, const double *l,const double *lpi,
                 const double *x,const double *xpi, const int m, double *xl);
void utmulx(double *y, const double *u, const double *xu, const int n);
void prpiutmulx(double *y, double *ypi, const double *u, const double *upi,
                const double *xu, const double *xupi, const int n);
void psdscaleK(double *y, const double *ud, const int *perm, const double *x,
            const coneK cK, const char transp, double *fwork);
void scaleK(double *y, double *dmult, const double *d, const double *ud,
            const double *qdetd, const int *perm, const double *x,
            const coneK cK, const char transp, const double *invdx,
            double *fwork);
#endif
