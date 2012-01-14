#include <math.h>

typedef struct {

  double      param;          /* 'c' parameter */
  double      hx, hy;        /* increment size in both directions */
  double      area;          /* area of the triangles */
} AppCtx;


int EPTorsLocalFunction(PetscInt coor[2], double x[4], double *f, void *ptr) {

  AppCtx *user = (AppCtx*)ptr;
  double c, hx, hy, area;

  double dvdx, dvdy, fquad, flin;

  c = user->param;
  hx = user->hx;
  hy = user->hy;
  area = user->area;

  /* 0 is 0,0; 1 is 1,0; 2 is 0,1; 3 is 1,1 */
  dvdx = (x[0] - x[1]) / hx;  /* lower triangle contribution */
  dvdy = (x[0] - x[2]) / hy;
  fquad = dvdx * dvdx + dvdy * dvdy;
  flin = x[0] + x[1] + x[2];
  dvdx = (x[3] - x[2]) / hx;  /* upper triangle contribution */
  dvdy = (x[3] - x[1]) / hy;
  fquad += dvdx * dvdx + dvdy * dvdy;
  flin += x[1] + x[2] + x[3];
  *f = area * (0.5 * fquad - (c * flin / 3.0)); /* 25 flops */

  return 0;
} /* LocalFunction */
