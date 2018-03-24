#if !defined(REFLECT_H)
#define REFLECT_H

void qtxq(double *x, const double *beta, const double *c,
          const int m, double *fwork);
void prpiqtxq(double *x, double *xpi, const double *beta, const double *c,
              const double *cpi, const int m, double *fwork);
void qxqt(double *x, const double *beta, const double *c,
          const int m, double *fwork);
void prpiqxqt(double *x, double *xpi, const double *beta, const double *c,
              const double *cpi, const int m, double *fwork);
#endif

