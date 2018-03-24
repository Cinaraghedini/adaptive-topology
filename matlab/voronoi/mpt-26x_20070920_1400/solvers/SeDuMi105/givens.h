#if !defined(GIVENS_H)
#define GIVENS_H

typedef struct{
  double x,y;} twodouble;

typedef struct{
  double x,xim,y;} tridouble;

void givensrot(double *z, const twodouble *g, const int n);
void givensrotuj(double *z, const twodouble *g, const int n);
void prpigivensrot(double *z, double *zpi, const tridouble *g, const int n);
void prpigivensrotuj(double *z, double *zpi, const tridouble *g, const int n);
#endif
