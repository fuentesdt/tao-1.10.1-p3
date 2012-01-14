#include <math.h>

typedef struct {
  InactiveDouble  param;          /* nonlinearity parameter */
  InactiveDouble  hx, hy, area;   /* increments and area of the triangle */
} AppCtx;

int CombLocalFunction(PetscInt coor[2], double x[4], double *f, void *ptr) {

  AppCtx *user = (AppCtx*)ptr;
  InactiveDouble lambda, hx, hy, area;

  double fquad, fexp, dvdx, dvdy;

  lambda = user->param;
  hx = user->hx;
  hy = user->hy;
  area = user->area;

  /* 0 is 0,0; 1 is 1,0; 2 is 0,1; 3 is 1,1 */
  dvdx = (x[0] - x[1]) / hx;  /* lower triangle contribution */
  dvdy = (x[0] - x[2]) / hy;
  fquad = dvdx * dvdx + dvdy * dvdy;
  fexp = exp(x[0]) + exp(x[1]) + exp(x[2]);
  dvdx = (x[3] - x[2]) / hx;  /* upper triangle contribution */
  dvdy = (x[3] - x[1]) / hy;
  fquad += dvdx * dvdx + dvdy * dvdy;
  fexp += exp(x[1]) + exp(x[2]) + exp(x[3]);
  *f = area * (0.5 * fquad - (lambda * fexp / 3.0)); /* 31 flops */

  return 0;
} /* LocalFunction */
