#include <math.h>

typedef struct {
  InactiveDouble      hx, hy;        /* increment size in both directions */
  InactiveDouble      area;          /* area of the triangles */
} AppCtx;

int MSurfLocalFunction(PetscInt coor[2], double x[4], double *f, void *ptr) {

  AppCtx *user = (AppCtx*)ptr;
  InactiveDouble hx, hy, area;

  double dvdx, dvdy, flow, fup;

  hx = user->hx;
  hy = user->hy;
  area = user->area;

  /* 0 is 0,0; 1 is 1,0; 2 is 0,1; 3 is 1,1 */
  dvdx = (x[0] - x[1]) / hx;  /* lower triangle contribution */
  dvdy = (x[0] - x[2]) / hy;
  flow = sqrt( 1 + dvdx * dvdx + dvdy * dvdy );

  dvdx = (x[3] - x[2]) / hx;  /* upper triangle contribution */
  dvdy = (x[3] - x[1]) / hy;
  fup = sqrt( 1 + dvdx * dvdx + dvdy * dvdy );

  *f = area * (flow + fup);  /* 18 flops */

  return 0;

} /* LocalFunction */
