#include <math.h>

/* typedef double InactiveDouble; */

#ifdef ad_GRAD_MAX
#undef ad_GRAD_MAX
#endif
#define ad_GRAD_MAX 16


typedef struct {
  int                 mx,my;
  InactiveDouble      hx, hy;        /* increment size in both directions */
  InactiveDouble      twopivornum; 
  InactiveDouble      tkappa;
} AppCtx;

int GinzLandLocalFunction(int coor[2], double xx[16], double *f, void *ptr) {

  AppCtx *user = (AppCtx*)ptr;
  InactiveDouble hx=user->hx, hy=user->hy;
  InactiveDouble tkappa = user->tkappa;
  InactiveDouble twopivornum=user->twopivornum;
  int i, mx=user->mx, my=user->my, sqn=mx*my;
  double arg;
  double d1,d2;
  double x1,x2,xy;
  double delsq, fcond,fkin,ffield;
  double x[16];


  for (i=0;i<16;i++){
    x[i]=xx[i];
  }

  if (coor[0]==mx-1){
      arg = (twopivornum*coor[1])/(my);
      d1=x[4]*cos(arg) - x[5]*sin(arg);
      d2=x[4]*sin(arg) + x[5]*cos(arg);

      x[4] = d1;
      x[5] = d2;
      x[7] = x[7] + (twopivornum) / (my*hy);
  }

  /*  Compute the Condensation Energy Density*/
  delsq = x[0]*x[0] + x[1]*x[1];
  fcond = (- delsq + (delsq*delsq)/2.0 );

  /*  Compute the Kinetic Energy Density. */
  d1=hx*x[2];
  x1 = x[4] - x[0]*cos(d1) + x[1]*sin(d1);
  x2=  x[5] - x[1]*cos(d1) - x[0]*sin(d1);
  fkin = (x1*x1 + x2*x2)/(hx*hx);
  d2=hy*x[3];
  x1 = x[8] - x[0]*cos(d2) + x[1]*sin(d2);
  x2 = x[9] - x[1]*cos(d2) - x[0]*sin(d2);
  fkin = fkin + (x1*x1 + x2*x2)/(hy*hy);

  /*  Compute the Magnetic Field Energy Density. */
  xy= (x[7] - x[3])/hx - (x[10] - x[2])/hy;
  ffield = xy*xy*(tkappa*tkappa);

  *f = (fcond + fkin + ffield)/sqn;

  return 0;
} /* LocalFunction */
