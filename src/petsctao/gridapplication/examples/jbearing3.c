/*$Id: jbearing.c, v 1.1 2002/08/08 11:35 lopezca@mauddib.mcs.anl.gov $*/

/* Program usage: mpirun -np <proc> jbearing [-help] [all TAO options] */

/*
  Include "tao.h" so we can use TAO solvers.
  petscda.h for distributed array
  ad_deriv.h for AD gradient
*/

#include "petscda.h"
#include "tao.h"
#include "taodaapplication.h"

static char help[] ="Pressure distribution in a Journal Bearing. \n\
This example is based on the problem DPJB from the MINPACK-2 test suite.\n\
This pressure journal bearing problem is an example of elliptic variational\n\
problem defined over a two dimensional rectangle. By discretizing the domain \n\
into triangular elements, the pressure surrounding the journal bearing is\n\
defined as the minimum of a quadratic function whose variables are bounded\n\
below by zero. The command line options are:\n\
  -ecc <ecc>, where <ecc> = epsilon parameter\n\
  -b <b>, where <b> = half the upper limit in the 2nd coordinate direction\n\
  -mx <xg>, where <xg> = number of grid points in the 1st coordinate direction\n\
  -my <yg>, where <yg> = number of grid points in the 2nd coordinate direction\n\
  -nlevels <nlevels>, where <nlevels> = number of levels in multigrid\n\
  -byelement, if computation is made by functions on rectangular elements\n\
  -adic, if AD is used (AD is not used by default)\n\n";

/*T
   Concepts: TAO - Solving a bounded minimization problem
   Routines: TaoInitialize(); TaoFinalize();
   Routines: TaoCreate(); TaoDestroy();
   Routines: DAApplicationCreate(); DAApplicationDestroy();
   Routines: DAAppSetVariableBoundsRoutine;
   Routines: DAAppSetElementObjectiveAndGradientRoutine();
   Routines: DAAppSetElementHessianRoutine();
   Routines: DAAppSetObjectiveAndGradientRoutine();
   Routines: DAAppSetADElementFunctionGradient();
   Routines: DAAppSetHessianRoutine();
   Routines: TaoSetOptions();
   Routines: TaoGetSolutionStatus(); TaoDAAppSolve();
   Routines: DAAppSetBeforeMonitor(); DAAppSetAfterMonitor
   Routines: DAAppGetSolution(); TaoView();
   Routines: DAAppGetInterpolationMatrix();
   Processors: n
T*/
 
/*
   User-defined application context - contains data needed by the
   application-provided call-back routines.
*/

int  ad_JBearLocalFunction(PetscInt[2] ,DERIV_TYPE[4], DERIV_TYPE *, void*);
typedef struct {

  InactiveDouble      *wq, *wl;      /* vectors with the parameters w_q(x) and w_l(x) */
  InactiveDouble      hx, hy;        /* increment size in both directions */
  InactiveDouble      area;          /* area of the triangles */

} ADFGCtx;


typedef struct {

  PetscReal      ecc;           /* epsilon value */
  PetscReal      b;             /* 0.5 * upper limit for 2nd variable */
  double      *wq, *wl;      /* vectors with the parameters w_q(x) and w_l(x) */
  double      hx, hy;        /* increment size in both directions */
  double      area;          /* area of the triangles */

  PetscInt    mx, my;        /* discretization including boundaries */

  ADFGCtx     fgctx;         /* Used only when an ADIC generated gradient is used */

} AppCtx;

/* User-defined routines found in this file */
static int AppCtxInitialize(void *ptr);
static int FormInitialGuess(DA, Vec);

static int JBearLocalFunctionGradient(PetscInt[2], double x[4], double *f, double g[4], void *ptr);
static int JBearLocalHessian(PetscInt[2], double x[4], double H[4][4], void *ptr);

static int WholeJBearFunctionGradient(TAO_APPLICATION,DA,Vec,double *,Vec,void*);
static int WholeJBearHessian(TAO_APPLICATION,DA,Vec,Mat,void*);

static int DASetBounds(TAO_APPLICATION, DA, Vec, Vec, void*);

static int MyGridMonitorBefore(TAO_APPLICATION, DA, PetscInt, void *);
static int MyGridMonitorAfter1(TAO_APPLICATION, DA, PetscInt, void *);

#undef __FUNCT__
#define __FUNCT__ "main"
int main( int argc, char **argv ) {

  PetscInt             info,iter;             /* used to check for functions returning nonzeros */
  PetscInt             nlevels;                                             /* multigrid levels */
  PetscInt             Nx,Ny;
  double          ff,gnorm;
  DA              DAarray[20];
  Vec             X;
  PetscTruth      flg, PreLoad = PETSC_TRUE;                                      /* flags */
  AppCtx          user;                                       /* user-defined work context */
  TaoMethod       method = "tao_gpcg";                              /* minimization method */
  TAO_SOLVER      tao;                                        /* TAO_SOLVER solver context */
  TAO_APPLICATION JBearApp;                                    /* The PETSc application */
  TaoTerminateReason reason;
  KSP ksp;
  PC  pc;

  /* Initialize TAO */
  PetscInitialize(&argc, &argv, (char *)0, help);
  TaoInitialize(&argc, &argv, (char *)0, help);

  PreLoadBegin(PreLoad,"Solve");
  
  info = AppCtxInitialize((void*)&user); CHKERRQ(info);
  
  nlevels=5;
  info = PetscOptionsGetInt(PETSC_NULL,"-nlevels",&nlevels,&flg); CHKERRQ(info);
  if (PreLoadIt == 0) {
    nlevels = 1; user.mx = 11; user.my = 11; }

  PetscPrintf(MPI_COMM_WORLD,"\n---- Journal Bearing Problem -----\n\n");

  /* Let PETSc determine the vector distribution */
  Nx = PETSC_DECIDE; Ny = PETSC_DECIDE;

  /* Create distributed array (DA) to manage parallel grid and vectors  */
  info = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,user.mx,
                    user.my,Nx,Ny,1,1,PETSC_NULL,PETSC_NULL,&DAarray[0]); CHKERRQ(info);
  for (iter=1;iter<nlevels;iter++){
    info = DARefine(DAarray[iter-1],PETSC_COMM_WORLD,&DAarray[iter]); CHKERRQ(info);
  }

  /* Create TAO solver and set desired solution method */
  info = TaoCreate(MPI_COMM_WORLD,method,&tao); CHKERRQ(info);
  info = TaoApplicationCreate(PETSC_COMM_WORLD, &JBearApp); CHKERRQ(info);
  info = TaoAppSetDAApp(JBearApp,DAarray,nlevels); CHKERRQ(info);

  /* Sets routine bounds evaluation */
  info = DAAppSetVariableBoundsRoutine(JBearApp,DASetBounds,(void *)&user); CHKERRQ(info);

  info = PetscOptionsHasName(TAO_NULL, "-byelement", &flg); CHKERRQ(info);
  if (flg) {

    /* Sets routines for function and gradient evaluation, element by element */
    info = PetscOptionsHasName(TAO_NULL, "-adic", &flg); CHKERRQ(info);
    if (flg) {
      info = DAAppSetADElementFunctionGradient(JBearApp,ad_JBearLocalFunction,248,(void *)&user.fgctx); CHKERRQ(info);
    } else {
      info = DAAppSetElementObjectiveAndGradientRoutine(JBearApp,JBearLocalFunctionGradient,63,(void *)&user); CHKERRQ(info);
    }
    /* Sets routines for Hessian evaluation, element by element */
    info = DAAppSetElementHessianRoutine(JBearApp,JBearLocalHessian,16,(void*)&user); CHKERRQ(info);

  } else {

    /* Sets routines for function and gradient evaluation, all in one routine */
    info = DAAppSetObjectiveAndGradientRoutine(JBearApp,WholeJBearFunctionGradient,(void *)&user); CHKERRQ(info);

    /* Sets routines for Hessian evaluation, all in one routine */
    info = DAAppSetHessianRoutine(JBearApp,WholeJBearHessian,(void*)&user); CHKERRQ(info);    
    
  }

  info = DAAppSetBeforeMonitor(JBearApp,MyGridMonitorBefore,(void*)&user); CHKERRQ(info);
  info = DAAppSetAfterMonitor(JBearApp,MyGridMonitorAfter1,(void*)&user); CHKERRQ(info);
  info = PetscOptionsHasName(TAO_NULL,"-tao_monitor", &flg); CHKERRQ(info);
  if (flg){
    info = DAAppPrintInterpolationError(JBearApp); CHKERRQ(info);
    info = DAAppPrintStageTimes(JBearApp); CHKERRQ(info);
  }
  /* Check for any tao command line options */
  info = TaoAppSetRelativeTolerance(JBearApp,1.0e-6); CHKERRQ(info);
  info = TaoSetTolerances(tao,0,0,0,0); CHKERRQ(info);
  info = TaoSetGradientTolerances(tao,0,0,0); CHKERRQ(info);

  info = TaoAppGetKSP(JBearApp,&ksp);CHKERRQ(info);
  info = KSPSetType(ksp,KSPCG); CHKERRQ(info);
  info = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,10);CHKERRQ(info);
  info = KSPGetPC(ksp,&pc);CHKERRQ(info);
  info = PCSetType(pc,PCBJACOBI);CHKERRQ(info);

  info = TaoSetOptions(JBearApp,tao); CHKERRQ(info);

  info = DAAppGetSolution(JBearApp,0,&X); CHKERRQ(info);
  info = FormInitialGuess(DAarray[0],X); CHKERRQ(info);
  info = DAAppSetInitialSolution(JBearApp,X); CHKERRQ(info);
  /* SOLVE THE APPLICATION */
  info = TaoDAAppSolve(JBearApp,tao);  CHKERRQ(info);

  /* Get information on termination */
  info = TaoGetSolutionStatus(tao,&iter,&ff,&gnorm,0,0,&reason); CHKERRQ(info);
  if (reason <= 0 ){
    PetscPrintf(MPI_COMM_WORLD,"Try a different TAO method, adjust some parameters, or check the function evaluation routines\n");
    PetscPrintf(MPI_COMM_WORLD," Iterations: %d,  Function Value: %4.2e, Residual: %4.2e \n",iter,ff,gnorm);
  }

  info = PetscOptionsHasName(PETSC_NULL,"-view_sol",&flg); CHKERRQ(info);
  if (flg){
    info = DAAppGetSolution(JBearApp,nlevels-1,&X); CHKERRQ(info);
    info=VecView(X,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(info);
  }

  /*  To View TAO solver information */
  // info = TaoView(tao); CHKERRQ(info);

  /* Free TAO data structures */
  info = TaoDestroy(tao); CHKERRQ(info);
  info = TaoAppDestroy(JBearApp); CHKERRQ(info);

  /* Free PETSc data structures */
  for (iter=0;iter<nlevels;iter++){
    info = DADestroy(DAarray[iter]); CHKERRQ(info);
  }

  PreLoadEnd();

  /* Finalize TAO */
  TaoFinalize();
  PetscFinalize();

  return 0;
} /* main */



/*----- The following two routines
   
  MyGridMonitorBefore    MyGridMonitorAfter

  help diplay info of iterations at every grid level -------*/

#undef __FUNCT__
#define __FUNCT__ "MyGridMonitorBefore"
static int MyGridMonitorBefore(TAO_APPLICATION myapp, DA da, PetscInt level, void *ctx) {

  AppCtx *user = (AppCtx*)ctx;
  int info;
  PetscInt mx,my;
  double t;

  info = DAGetInfo(da,PETSC_NULL,&mx,&my,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,
		   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(info);
  user->mx = mx;
  user->my = my;

  user->hx = (2.0 * 3.14159265358979) / (user->mx - 1);
  user->hy = (2.0 * user->b) / (user->my - 1);
  user->area = 0.5 * user->hx * user->hy;

  user->wq = new double[user->mx];
  user->wl = new double[user->mx];
  for (PetscInt i=0; i<user->mx; i++) {
    t = 1.0 + user->ecc * cos(i*user->hx);
    user->wq[i] = t*t*t;
    user->wl[i] = user->ecc * sin(i*user->hx);
  }
  info = PetscLogFlops(8 + 7 * user->mx); CHKERRQ(info);

  user->fgctx.hx   = user->hx;
  user->fgctx.hy   = user->hy;
  user->fgctx.area = user->area;
  user->fgctx.wq = user->wq;
  user->fgctx.wl = user->wl;

  PetscPrintf(MPI_COMM_WORLD,"Grid: %d,    mx: %d     my: %d   \n",level,mx,my);

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "MyGridMonitorAfter1"
static int MyGridMonitorAfter1(TAO_APPLICATION myapp, DA da, PetscInt level, void *ctx){
  
  AppCtx *user = (AppCtx*)ctx;
  delete [] user->wq; delete [] user->wl;
  return 0;
}



/*------- USER-DEFINED: initialize the application context information -------*/

#undef __FUNCT__
#define __FUNCT__ "AppCtxInitialize"
/*
  AppCtxInitialize - Sets initial values for the application context parameters

  Input:
    ptr - void user-defined application context

  Output:
    ptr - user-defined application context with the default or user-provided
             parameters
*/
static int AppCtxInitialize(void *ptr) {

  AppCtx *user = (AppCtx*)ptr;
  PetscTruth    flg;            /* flag for PETSc calls */
  int info;

  /* Specify default parameters */
  user->ecc = 0.1;
  user->b = 10.0;
  user->mx = user->my = 11;

  /* Check for command line arguments that override defaults */
  info = PetscOptionsGetReal(TAO_NULL, "-ecc", &user->ecc, &flg); CHKERRQ(info);
  info = PetscOptionsGetReal(TAO_NULL, "-b", &user->b, &flg); CHKERRQ(info);
  info = PetscOptionsGetInt(TAO_NULL, "-mx", &user->mx, &flg); CHKERRQ(info);
  info = PetscOptionsGetInt(TAO_NULL, "-my", &user->my, &flg); CHKERRQ(info);

  return 0;
} /* AppCtxInitialize */


#undef __FUNCT__
#define __FUNCT__ "FormInitialGuess"
static int FormInitialGuess(DA da, Vec X)
{
  int info;
  PetscInt    i, j, mx;
  PetscInt    xs, ys, xm, ym, xe, ye;
  PetscReal hx, val;
  double **x;

  /* Get local mesh boundaries */
  info = DAGetInfo(da,PETSC_NULL,&mx,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(info);
  hx = 2.0*4.0*atan(1.0)/(mx-1);

  info = DAGetCorners(da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL); CHKERRQ(info);
  xe = xs+xm; ye = ys+ym;

  info = DAVecGetArray(da, X, (void**)&x); CHKERRQ(info);
  /* Compute initial guess over locally owned part of mesh */
  for (j=ys; j<ye; j++) {  /*  for (j=0; j<my; j++) */
    for (i=xs; i<xe; i++) {  /*  for (i=0; i<mx; i++) */
      val = PetscMax(sin((i+1.0)*hx),0.0);
      x[j][i] = val;
      x[j][i] = 0;
    }
  }
  info = DAVecRestoreArray(da, X, (void**)&x); CHKERRQ(info);

  return 0;
}


/*------- USER-DEFINED: set the upper and lower bounds for the variables  -------*/
#undef __FUNCT__
#define __FUNCT__ "DASetBounds"
/*
  FormBounds - Forms bounds on the variables

  Input:
    user - user-defined application context

  Output:
    XL - vector of lower bounds
    XU - vector of upper bounds

*/
static int DASetBounds(TAO_APPLICATION daapplication, DA da, Vec XL, Vec XU, void *ptr) {

  AppCtx *user = (AppCtx*)ptr;
  int info;
  PetscInt i, j, mx, my;
  PetscInt xs, xm, ys, ym;
  double **xl, **xu;

  mx = user->mx;
  my = user->my;

  info = DAVecGetArray(da, XL, (void**)&xl); CHKERRQ(info);
  info = DAVecGetArray(da, XU, (void**)&xu); CHKERRQ(info);
  info = DAGetCorners(da, &xs, &ys, TAO_NULL, &xm, &ym, TAO_NULL); CHKERRQ(info);

  for (j = ys; j < ys+ym; j++){
    for (i = xs; i < xs+xm; i++){
      xl[j][i] = 0.0;
      if (i == 0 || j == 0 || i == mx - 1 || j == my - 1) {
        xu[j][i] = 0.0;
      } else {
        xu[j][i] = TAO_INFINITY;
      }
    }
  }

  info = DAVecRestoreArray(da, XL, (void**)&xl); CHKERRQ(info);
  info = DAVecRestoreArray(da, XU, (void**)&xu); CHKERRQ(info);

  return 0;

} /* DASetBounds */



#undef __FUNCT__
#define __FUNCT__ "JBearLocalFunctionGradient"
/*
  JBearLocalFunctionGradient - Evaluates function and gradient over the 
      local rectangular element

  Input:
    coor - vector with the indices of the position of current element
             in the first, second and third directions
    x - current point (values over the current rectangular element)
    ptr - user-defined application context

  Output:
    f - value of the objective funtion at the local rectangular element
    g - gradient of the local function

*/
static int JBearLocalFunctionGradient(PetscInt coor[2], double x[4], double *f, double g[4], void *ptr) {

  AppCtx *user = (AppCtx*)ptr;

  double avgWq, sqGrad, avgWV, fl, fu;
  double hx, hy, area, aread3, *wq, *wl;
  double dvdx, dvdy;
  PetscInt i;

  hx = user->hx;
  hy = user->hy;
  area = user->area;
  aread3 = area / 3.0;
  wq = user->wq;
  wl = user->wl;
  i = coor[0];

  /* lower triangle contribution */
  dvdx = (x[0] - x[1]) / hx;
  dvdy = (x[0] - x[2]) / hy;
  sqGrad = dvdx * dvdx + dvdy * dvdy;
  avgWq = (2.0 * wq[i] + wq[i+1]) / 6.0;
  avgWV = (wl[i]*x[0] + wl[i+1]*x[1] + wl[i]*x[2]) / 3.0;
  fl = avgWq * sqGrad - avgWV;

  dvdx = dvdx * hy * avgWq;
  dvdy = dvdy * hx * avgWq;
  g[0] = ( dvdx + dvdy ) - wl[i] * aread3;
  g[1] = ( -dvdx ) - wl[i+1] * aread3;
  g[2] = ( -dvdy ) - wl[i] * aread3;

  /* upper triangle contribution */
  dvdx = (x[3] - x[2]) / hx; 
  dvdy = (x[3] - x[1]) / hy;
  sqGrad = dvdx * dvdx + dvdy * dvdy;
  avgWq = (2.0 * wq[i+1] + wq[i]) / 6.0;
  avgWV = (wl[i+1]*x[1] + wl[i]*x[2] + wl[i+1]*x[3]) / 3.0;
  fu = avgWq * sqGrad - avgWV;

  dvdx = dvdx * hy * avgWq;
  dvdy = dvdy * hx * avgWq;
  g[1] += (-dvdy) - wl[i+1] * aread3;
  g[2] +=  (-dvdx) - wl[i] * aread3;
  g[3] = ( dvdx + dvdy ) - wl[i+1] * aread3;

  *f = area * (fl + fu);

  return 0;
} /* JBearLocalFunctionGradient */


/*------- USER-DEFINED: routine to evaluate the Hessian
           at a local (rectangular element) level       -------*/
#undef __FUNCT__
#define __FUNCT__ "JBearLocalHessian"
/*
  JBearLocalHessian - Computes the Hessian of the local (partial) function
         defined over the current rectangle

  Input:
    coor - vector with the indices of the position of current element
             in the first, second and third directions
    x - current local solution (over the rectangle only)
    ptr - user-defined application context

  Output:
    H - Hessian matrix of the local function (wrt the four
           points of the rectangle only)

*/
static int JBearLocalHessian(PetscInt coor[2], double x[4], double H[4][4],void *ptr) {

  AppCtx *user = (AppCtx*)ptr;
  double wql, wqu;
  double hx, hy, dydx, dxdy;
  double *wq;
  double wqldydx,wqldxdy,wqudydx,wqudxdy;
  PetscInt i;

  hx = user->hx;
  hy = user->hy;
  wq = user->wq;
  i = coor[0];

  dxdy = hx / hy;
  dydx = hy / hx;
  wql = (2.0 * wq[i] + wq[i+1]) / 6.0;
  wqu = (wq[i] + 2.0 * wq[i+1]) / 6.0;
  wqldydx = wql * dydx;
  wqldxdy = wql * dxdy;
  wqudydx = wqu * dydx;
  wqudxdy = wqu * dxdy;

          /* Hessian contribution at 0,0 */
  H[0][0] = wqldxdy + wqldydx;
  H[0][1] = H[1][0] = -wqldydx;
  H[0][2] = H[2][0] = -wqldxdy;
  H[0][3] = H[3][0] = 0.0;

          /* Hessian contribution at 1,0 */
  H[1][1] = wqldydx + wqudxdy;
  H[1][2] = H[2][1] = 0.0;
  H[1][3] = H[3][1] = -wqudxdy; 

          /* Hessian contribution at 0,1 */
  H[2][2] = wqldxdy + wqudydx;
  H[2][3] = H[3][2] = -wqudydx; 

          /* Hessian contribution at 1,1 */
  H[3][3] = wqudydx + wqudxdy;

  return 0;

} /* JBearLocalHessian */


/*------- USER-DEFINED: routine to evaluate the function 
          and gradient at the whole grid             -------*/

#undef __FUNCT__
#define __FUNCT__ "WholeJBearFunctionGradient"
/*
  WholeJBearFunctionGradient - Evaluates function and gradient over the 
      whole grid

  Input:
    daapplication - TAO application object
    da  - distributed array
    X   - the current point, at which the function and gradient are evaluated
    ptr - user-defined application context

  Output:
    f - value of the objective funtion at X
    G - gradient at X
*/
static int WholeJBearFunctionGradient(TAO_APPLICATION daapplication, DA da, Vec X, double *f, Vec G, void *ptr) {

  AppCtx *user = (AppCtx*)ptr;
  Vec localX, localG;
  int info;
  PetscInt i, j;
  PetscInt xs, xm, gxs, gxm, xe, ys, ym, gys, gym, ye;
  double **x, **g;
  double floc = 0.0;
  PetscScalar zero = 0.0;

  double avgWq, sqGrad, avgWV, fl, fu;
  double hx, hy, area, aread3, *wq, *wl;
  double dvdx, dvdy;

  hx = user->hx;
  hy = user->hy;
  area = user->area;
  aread3= area/3.0;
  wq = user->wq;
  wl = user->wl;

  info = DAGetLocalVector(da, &localX); CHKERRQ(info);
  info = DAGetLocalVector(da, &localG); CHKERRQ(info);
  info = VecSet(G, zero); CHKERRQ(info);
  info = VecSet(localG, zero); CHKERRQ(info);

  info = DAGlobalToLocalBegin(da, X, INSERT_VALUES, localX); CHKERRQ(info);
  info = DAGlobalToLocalEnd(da, X, INSERT_VALUES, localX); CHKERRQ(info);

  info = DAVecGetArray(da, localX, (void**)&x); CHKERRQ(info);
  info = DAVecGetArray(da, localG, (void**)&g); CHKERRQ(info);

  info = DAGetCorners(da, &xs, &ys, TAO_NULL, &xm, &ym, TAO_NULL); CHKERRQ(info);
  info = DAGetGhostCorners(da, &gxs, &gys, TAO_NULL, &gxm, &gym, TAO_NULL); CHKERRQ(info);

  xe = gxs + gxm - 1;
  ye = gys + gym - 1;
  for (j = ys; j < ye; j++) {
    for (i = xs; i < xe; i++) {

      /* lower triangle contribution */
      dvdx = (x[j][i] - x[j][i+1]) / hx;
      dvdy = (x[j][i] - x[j+1][i]) / hy;
      sqGrad = dvdx * dvdx + dvdy * dvdy;
      avgWq = (2.0 * wq[i] + wq[i+1]) / 6.0;
      avgWV = (wl[i]*x[j][i] + wl[i+1]*x[j][i+1] + wl[i]*x[j+1][i]) / 3.0;
      fl = avgWq * sqGrad - avgWV;

      dvdx = dvdx * hy * avgWq;
      dvdy = dvdy * hx * avgWq;
      g[j][i] +=  ( dvdx + dvdy ) - wl[i] * aread3;
      g[j][i+1] += ( -dvdx ) - wl[i+1] * aread3;
      g[j+1][i] += ( -dvdy ) - wl[i] * aread3;

      /* upper triangle contribution */
      dvdx = (x[j+1][i+1] - x[j+1][i]) / hx;
      dvdy = (x[j+1][i+1] - x[j][i+1]) / hy;
      sqGrad = dvdx * dvdx + dvdy * dvdy;
      avgWq = (2.0 * wq[i+1] + wq[i]) / 6.0;
      avgWV = (wl[i+1]*x[j][i+1] + wl[i]*x[j+1][i] + wl[i+1]*x[j+1][i+1]) / 3.0;
      fu = avgWq * sqGrad - avgWV;

      dvdx = dvdx * hy * avgWq;
      dvdy = dvdy * hx * avgWq;
      g[j][i+1] += (-dvdy) - wl[i+1] * aread3;
      g[j+1][i] +=  (-dvdx) - wl[i] * aread3;
      g[j+1][i+1] += ( dvdx + dvdy ) - wl[i+1] * aread3;

      floc += area * (fl + fu);

    }
  }

  info = MPI_Allreduce(&floc, f, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); CHKERRQ(info);

  info = DAVecRestoreArray(da, localX, (void**)&x); CHKERRQ(info);
  info = DAVecRestoreArray(da, localG, (void**)&g); CHKERRQ(info);

  info = DALocalToGlobalBegin(da, localG, G); CHKERRQ(info);
  info = DALocalToGlobalEnd(da, localG, G); CHKERRQ(info);

  info = DARestoreLocalVector(da, &localX); CHKERRQ(info);
  info = DARestoreLocalVector(da, &localG); CHKERRQ(info);

  info = PetscLogFlops((xe-xs) * (ye-ys) * 67 + 1); CHKERRQ(info);
  return 0;
} /* WholeJBearFunctionGradient */


#undef __FUNCT__
#define __FUNCT__ "WholeJBearHessian"
/*
  WholeJBearHessian - Evaluates Hessian over the whole grid

  Input:
    daapplication - TAO application object
    da  - distributed array
    X   - the current point, at which the function and gradient are evaluated
    ptr - user-defined application context

  Output:
    H - Hessian at X
*/
static int WholeJBearHessian(TAO_APPLICATION daapplication, DA da, Vec X, Mat H, void *ptr) {

  AppCtx *user = (AppCtx*)ptr;
  int info;
  PetscInt i, j, ind[4];
  PetscInt xs, xm, gxs, gxm, xe, ys, ym, gys, gym, ye;
  double smallH[4][4];
  double wql, wqu;
  double dydx, dxdy;
  double *wq;
  double wqldydx,wqldxdy,wqudydx,wqudxdy;
  PetscTruth assembled;

  wq = user->wq;
  dydx = user->hy / user->hx;
  dxdy = user->hx / user->hy;

  info = MatAssembled(H,&assembled); CHKERRQ(info);
  if (assembled){info = MatZeroEntries(H);  CHKERRQ(info);}


  info = DAGetCorners(da, &xs, &ys, TAO_NULL, &xm, &ym, TAO_NULL); CHKERRQ(info);
  info = DAGetGhostCorners(da, &gxs, &gys, TAO_NULL, &gxm, &gym, TAO_NULL); CHKERRQ(info);

  xe = gxs + gxm - 1;
  ye = gys + gym - 1;
  for (j = ys; j < ye; j++) {
    for (i = xs; i < xe; i++) {

      wql = (2.0 * wq[i] + wq[i+1]) / 6.0;
      wqu = (wq[i] + 2.0 * wq[i+1]) / 6.0;

      wqldydx = wql * dydx;
      wqldxdy = wql * dxdy;
      wqudydx = wqu * dydx;
      wqudxdy = wqu * dxdy;

          /* Hessian contribution at 0,0 */
      smallH[0][0] = wqldxdy + wqldydx;
      smallH[0][1] = smallH[1][0] = -wqldydx;
      smallH[0][2] = smallH[2][0] = -wqldxdy;
      smallH[0][3] = smallH[3][0] = 0.0;

          /* Hessian contribution at 1,0 */
      smallH[1][1] = (wqldydx + wqudxdy);
      smallH[1][2] = smallH[2][1] = 0.0;
      smallH[1][3] = smallH[3][1] = -wqudxdy;

          /* Hessian contribution at 0,1 */
      smallH[2][2] = (wqldxdy + wqudydx);
      smallH[2][3] = smallH[3][2] = -wqudydx;

          /* Hessian contribution at 1,1 */
      smallH[3][3] = wqudydx + wqudxdy;

      ind[0] = (j-gys) * gxm + (i-gxs);
      ind[1] = ind[0] + 1;
      ind[2] = ind[0] + gxm;
      ind[3] = ind[2] + 1;
      info = MatSetValuesLocal(H,4,ind,4,ind,(PetscScalar*)smallH,ADD_VALUES); CHKERRQ(info);

    }
  }


  info = MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  info = MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  info = MatSetOption(H, MAT_SYMMETRIC, PETSC_TRUE); CHKERRQ(info);


  info = PetscLogFlops((xe-xs) * (ye-ys) * 14 + 2); CHKERRQ(info);
  return 0;

} /* WholeJBearHessian */
