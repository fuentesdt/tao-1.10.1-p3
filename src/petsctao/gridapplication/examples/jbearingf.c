#include <math.h>

typedef struct {

  InactiveDouble      *wq, *wl;      /* vectors with the parameters w_q(x) and w_l(x) */
  InactiveDouble      hx, hy;        /* increment size in both directions */
  InactiveDouble      area;          /* area of the triangles */

} AppCtx;


int JBearLocalFunction(PetscInt coor[2], double x[4], double *f, void *ptr) {

  AppCtx *user = (AppCtx*)ptr;
  InactiveDouble hx, hy, area;
  InactiveDouble *wq, *wl;

  double avgWq, sqGrad, avgWV;
  double dvdx, dvdy, fl, fu;
  PetscInt i;

  hx = user->hx;
  hy = user->hy;
  area = user->area;
  wq = user->wq;
  wl = user->wl;
  i = coor[0];

  /* 0 is 0,0; 1 is 1,0; 2 is 0,1; 3 is 1,1 */
  dvdx = (x[0] - x[1]) / hx;  /* lower triangle contribution */
  dvdy = (x[0] - x[2]) / hy;
  sqGrad = dvdx * dvdx + dvdy * dvdy;
  avgWq = (2.0 * wq[i] + wq[i+1]) / 3.0;
  avgWV = (wl[i]*x[0] + wl[i+1]*x[1] + wl[i]*x[2]) / 3.0;
  fl = 0.5 * avgWq * sqGrad - avgWV;

  dvdx = (x[3] - x[2]) / hx;  /* upper triangle contribution */
  dvdy = (x[3] - x[1]) / hy;
  sqGrad = dvdx * dvdx + dvdy * dvdy;
  avgWq = (2.0 * wq[i+1] + wq[i]) / 3.0;
  avgWV = (wl[i+1]*x[1] + wl[i]*x[2] + wl[i+1]*x[3]) / 3.0;
  fu = 0.5 * avgWq * sqGrad - avgWV;

  *f = area * (fl + fu);   /* 40 flops */

  return 0;
} /* LocalFunction */
