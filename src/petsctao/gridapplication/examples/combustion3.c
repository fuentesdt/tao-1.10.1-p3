/*$Id: combustion.c,v 1.1 2002/08/08 9:39 lopezca@mauddib.mcs.anl.gov $*/

/* Program usage: mpirun -np <proc> combustion [-help] [all TAO options] */

/*
  Include "tao.h" so we can use TAO solvers.
  petscda.h for distributed array
  ad_deriv.h for AD gradient
*/

#include "petscda.h"
#include "tao.h"
#include "taodaapplication.h"

static char help[] = "Steady-State Combustion.\n\
We solve the Steady-State Combustion problem (MINPACK-2 test suite) in a 2D\n\
rectangular domain, using distributed arrays (DAs) to partition the parallel grid.\n\
The command line options include:\n\
  -mx <xg>, where <xg> = number of grid points in the 1st coordinate direction\n\
  -my <yg>, where <yg> = number of grid points in the 2nd coordinate direction\n\
  -nlevels <nlevels>, where <nlevels> = number of levels in multigrid\n\
  -byelement, if computation is made by functions on rectangular elements\n\
  -adic, if AD is used (AD is not used by default)\n\
  -par <parameter>, where <parameter> indicates the problem's nonlinearity\n\
     parameter lambda (0 <= par <= 6.81)\n\n";

/*T
   Concepts: TAO - Solving a bounded minimization problem
   Routines: TaoInitialize(); TaoFinalize();
   Routines: TaoCreate(); TaoDestroy();
   Routines: DAApplicationCreate(); DAApplicationDestroy();
   Routines: DAAppSetVariableBoundsRoutine();
   Routines: DAAppSetElementObjectiveAndGradientRoutine();
   Routines: DAAppSetElementHessianRoutine();
   Routines: DAAppSetObjectiveAndGradientRoutine();
   Routines: DAAppSetADElementFunctionGradient();
   Routines: DAAppSetHessianRoutine();
   Routines: TaoAppSetOptions();
   Routines: TaoGetSolutionStatus(); TaoDAAppSolve();
   Routines: DAAppSetMonitor(); TaoView();
   Routines: DAAppGetSolution();
   Routines: DAAppGetInterpolationMatrix();
   Processors: n
T*/

/*
   User-defined application context - contains data needed by the
   application-provided call-back routines.
*/

typedef struct {
  InactiveDouble      param;
  InactiveDouble      hx, hy;        /* increment size in both directions */
  InactiveDouble      area;          /* area of the triangles */
} ADFGCtx;

typedef struct {
  PetscReal  param;          /* nonlinearity parameter */
  double  hx, hy, area;   /* increments and area of the triangle */
  PetscInt     mx, my;         /* discretization including boundaries */
  ADFGCtx fgctx;          /* Used only when an ADIC generated gradient is used */
} AppCtx;
int ad_CombLocalFunction(PetscInt[2], DERIV_TYPE[4], DERIV_TYPE*, void*);

/* User-defined routines foun in this file */
static int AppCtxInitialize(void *ptr);
static int FormInitialGuess(DA, Vec, AppCtx*);

static int CombLocalFunctionGradient(PetscInt[3], double x[4], double *f, double g[4], void *ptr);
static int WholeCombFunctionGradient(TAO_APPLICATION,DA,Vec,double *,Vec,void*);

static int CombLocalHessian(PetscInt[3], double x[4], double H[4][4], void *ptr);
static int WholeCombHessian(TAO_APPLICATION,DA,Vec,Mat,void*);

static int DAFixBoundary(TAO_APPLICATION, DA, Vec, Vec, void*);

static int MyGridMonitorBefore(TAO_APPLICATION, DA, PetscInt, void *);

#undef __FUNCT__
#define __FUNCT__ "main"

int main( int argc, char **argv ) {

  int             info;           /* used to check for functions returning nonzeros */
  PetscInt        Nx,Ny,iter;
  PetscInt        nlevels;                                           /* multigrid levels */
  double          ff,gnorm;
  DA              DAarray[20];
  Vec             X;
  KSP             ksp;
  PetscTruth      flg, PreLoad = PETSC_TRUE;                                    /* flags */
  AppCtx          user;                                     /* user-defined work context */
  TaoMethod       method = "tao_tron";                            /* minimization method */
  TAO_SOLVER      tao;                                      /* TAO_SOLVER solver context */
  TAO_APPLICATION CombApp;                                   /* The PETSc application */
  TaoTerminateReason reason;

  /* Initialize TAO */
  PetscInitialize(&argc, &argv, (char *)0, help);
  TaoInitialize(&argc, &argv, (char *)0, help);

  PreLoadBegin(PreLoad,"Solve");
  
  info = AppCtxInitialize((void*)&user); CHKERRQ(info);

  nlevels=5;
  info = PetscOptionsGetInt(PETSC_NULL,"-nlevels",&nlevels,&flg); CHKERRQ(info);
  if (PreLoadIt == 0) {
    nlevels = 1; user.mx = 11; user.my = 11;}

  PetscPrintf(MPI_COMM_WORLD,"\n---- Steady-State Combustion Problem -----\n\n");

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

  info = TaoApplicationCreate(PETSC_COMM_WORLD, &CombApp); CHKERRQ(info);
  info = TaoAppSetDAApp(CombApp,DAarray,nlevels); CHKERRQ(info);

  /* Sets routines for function, gradient and bounds evaluation */
  info = DAAppSetVariableBoundsRoutine(CombApp,DAFixBoundary,(void *)&user); CHKERRQ(info);

  info = PetscOptionsHasName(TAO_NULL, "-byelement", &flg); CHKERRQ(info);
  if (flg) {

    /* Sets routines for function and gradient evaluation, element by element */
    info = PetscOptionsHasName(TAO_NULL, "-adic", &flg); CHKERRQ(info);
    if (flg) {
      info = DAAppSetADElementFunctionGradient(CombApp,ad_CombLocalFunction,228,(void *)&user.fgctx); CHKERRQ(info);
    } else {
      info = DAAppSetElementObjectiveAndGradientRoutine(CombApp,CombLocalFunctionGradient,51,(void *)&user); CHKERRQ(info);
    }
    /* Sets routines for Hessian evaluation, element by element */
    info = DAAppSetElementHessianRoutine(CombApp,CombLocalHessian,21,(void*)&user); CHKERRQ(info);

  } else {

    /* Sets routines for function and gradient evaluation, all in one routine */
    info = DAAppSetObjectiveAndGradientRoutine(CombApp,WholeCombFunctionGradient,(void *)&user); CHKERRQ(info);

    /* Sets routines for Hessian evaluation, all in one routine */
    info = DAAppSetHessianRoutine(CombApp,WholeCombHessian,(void*)&user); CHKERRQ(info);    
    
  }

  info = DAAppSetBeforeMonitor(CombApp,MyGridMonitorBefore,(void*)&user); CHKERRQ(info);
  info = PetscOptionsHasName(TAO_NULL,"-tao_monitor", &flg); CHKERRQ(info);
  if (flg){
    info = DAAppPrintInterpolationError(CombApp); CHKERRQ(info);
    info = DAAppPrintStageTimes(CombApp); CHKERRQ(info);
  }


  info = TaoAppSetRelativeTolerance(CombApp,1.0e-6); CHKERRQ(info);
  info = TaoSetTolerances(tao,0,0,0,0); CHKERRQ(info);
  info = TaoSetGradientTolerances(tao,0,0,0); CHKERRQ(info);

  info = TaoAppGetKSP(CombApp,&ksp);CHKERRQ(info);
  info = KSPSetType(ksp,KSPCG); CHKERRQ(info);
  info = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,100);CHKERRQ(info);

  /* Check for any tao command line options */
  info = TaoSetOptions(CombApp, tao); CHKERRQ(info);

  info = DAAppGetSolution(CombApp,0,&X); CHKERRQ(info);
  info = FormInitialGuess(DAarray[0],X,&user); CHKERRQ(info);
  info = DAAppSetInitialSolution(CombApp,X); CHKERRQ(info);

  /* SOLVE THE APPLICATION */
  info = TaoDAAppSolve(CombApp, tao);  CHKERRQ(info);

  /* Get information on termination */
  info = TaoGetSolutionStatus(tao,&iter,&ff,&gnorm,0,0,&reason); CHKERRQ(info);
  if (reason <= 0 ){
    PetscPrintf(MPI_COMM_WORLD,"Try a different TAO method, adjust some parameters, or check the function evaluation routines\n");
    PetscPrintf(MPI_COMM_WORLD," Iterations: %d,  Function Value: %4.2e, Residual: %4.2e \n",iter,ff,gnorm);
  }

  info = PetscOptionsHasName(PETSC_NULL,"-view_sol",&flg); CHKERRQ(info);
  if (flg){
    info = DAAppGetSolution(CombApp,nlevels-1,&X); CHKERRQ(info);
    info=VecView(X,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(info);
  }

  /*  To View TAO solver information */
  // info = TaoView(tao); CHKERRQ(info);

  /* Free TAO data structures */
  info = TaoDestroy(tao); CHKERRQ(info);
  info = TaoAppDestroy(CombApp); CHKERRQ(info);

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

  info = DAGetInfo(da,PETSC_NULL,&mx,&my,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,
		   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(info);
  user->mx = mx;
  user->my = my;
  user->hx = 1.0 / (user->mx - 1);
  user->hy = 1.0 / (user->my - 1);
  user->area = 0.5 * user->hx * user->hy;
  user->fgctx.hx   = user->hx;
  user->fgctx.hy   = user->hy;
  user->fgctx.area = user->area;
  user->fgctx.param = user->param;

  PetscPrintf(MPI_COMM_WORLD,"Grid: %d,    mx: %d     my: %d   \n",level,mx,my);
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
  PetscReal     LambdaMax = 6.81, LambdaMin = 0.0;  /* bounds on parameter lambda */
  PetscTruth    flg;            /* flag for PETSc calls */
  int info;

  /* Specify dimension of the problem */
  user->param = 5.0;
  user->mx = 11;
  user->my = 11;

  /* Check for any command line arguments that override defaults */
  info = PetscOptionsGetReal(TAO_NULL, "-par", &user->param, &flg); CHKERRQ(info);
  if (user->param >= LambdaMax || user->param <= LambdaMin) {
    SETERRQ(1,"Lambda is out of range.");
  }
  info = PetscOptionsGetInt(PETSC_NULL,"-mx",&user->mx,&flg); CHKERRQ(info);
  info = PetscOptionsGetInt(PETSC_NULL,"-my",&user->my,&flg); CHKERRQ(info);

  user->hx = 1.0 / (user->mx - 1);
  user->hy = 1.0 / (user->my - 1);
  user->area = 0.5 * user->hx * user->hy;
  info = PetscLogFlops(6); CHKERRQ(info);

  return 0;
} /* AppCtxInitialize */



#undef __FUNCT__
#define __FUNCT__ "FormInitialGuess"
static int FormInitialGuess(DA da, Vec X, AppCtx *ctx)
{
  int    info;
  PetscInt i, j, mx, my;
  PetscInt xs, ys, xm, ym, xe, ye;
  PetscReal hx, hy, temp, val, lambda;
  double **x;

  lambda = ctx->param;
  lambda = lambda/(lambda+1.0);

  /* Get local mesh boundaries */
  info = DAGetInfo(da,PETSC_NULL,&mx,&my,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(info);
  hx = 1.0/(mx-1);  hy = 1.0/(my-1);

  info = DAGetCorners(da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL); CHKERRQ(info);
  xe = xs+xm; ye = ys+ym;

  info = DAVecGetArray(da, X, (void**)&x); CHKERRQ(info);
  /* Compute initial guess over locally owned part of mesh */
  for (j=ys; j<ye; j++) {  /*  for (j=0; j<my; j++) */
    temp = PetscMin(j+1,my-j)*hy;
    for (i=xs; i<xe; i++) {  /*  for (i=0; i<mx; i++) */
      val = lambda*sqrt(PetscMin((PetscMin(i+1,mx-i))*hx,temp));
      x[j][i] = val;
    }
  }
  info = DAVecRestoreArray(da, X, (void**)&x); CHKERRQ(info);

  return 0;
}


/*------- USER-DEFINED: set the upper and lower bounds for the variables  -------*/

#undef __FUNCT__
#define __FUNCT__ "DAFixBoundary"
/*
  FormBounds - Forms bounds on the variables

  Input:
    user - user-defined application context

  Output:
    XL - vector of lower bounds
    XU - vector of upper bounds

*/
static int DAFixBoundary(TAO_APPLICATION daapplication, DA da, Vec XL, Vec XU, void *ptr)
{
  AppCtx *user = (AppCtx*)ptr;
  int info;
  PetscInt i, j, mx, my, xs, xm, ys, ym;
  double lb = -TAO_INFINITY;
  double ub = TAO_INFINITY;
  double **xl, **xu;

  mx = user->mx;  /* number of points including the boundary */
  my = user->my;

  info = DAVecGetArray(da, XL, (void**)&xl); CHKERRQ(info);
  info = DAVecGetArray(da, XU, (void**)&xu); CHKERRQ(info);
  info = DAGetCorners(da, &xs, &ys, TAO_NULL, &xm, &ym, TAO_NULL); CHKERRQ(info);

  /* Compute initial guess over locally owned part of the grid */
  for (j = ys; j < ys+ym; j++){
    for (i = xs; i < xs+xm; i++){
      if (i == 0 || j == 0 || i == mx - 1 || j == my - 1) {
        xl[j][i] = xu[j][i] = 0.0;
      } else {
        xl[j][i] = lb;
        xu[j][i] = ub;
      }
    }
  }

  info = DAVecRestoreArray(da, XL, (void**)&xl); CHKERRQ(info);
  info = DAVecRestoreArray(da, XU, (void**)&xu); CHKERRQ(info);

  return 0;
} /* DAFixBoundary */


/*------- USER-DEFINED: routine to evaluate the function and gradient
           at a local (rectangular element) level              -------*/

#undef __FUNCT__
#define __FUNCT__ "CombLocalFunctionGradient"
/*
  CombLocalFunctionGradient - Evaluates function and gradient over the 
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
static int CombLocalFunctionGradient(PetscInt coor[2], double x[4], double *f, double g[4], void *ptr) {

  AppCtx *user = (AppCtx*)ptr;

  double lambdad3, hx, hy, area;
  double fquad, fexp, dvdx, dvdy;

  lambdad3 = user->param / 3.0;
  hx = user->hx;
  hy = user->hy;
  area = user->area;

  /* lower triangle contribution */
  dvdx = (x[0] - x[1]) / hx;
  dvdy = (x[0] - x[2]) / hy;
  fquad = dvdx * dvdx + dvdy * dvdy;
  fexp = exp(x[0]) + exp(x[1]) + exp(x[2]);

  dvdx = 0.5 * dvdx * hy;
  dvdy = 0.5 * dvdy * hx;
  g[0] = dvdx + dvdy - exp(x[0]) * area * lambdad3;
  g[1] = -dvdx - 2.0 * exp(x[1]) * area * lambdad3;
  g[2] = -dvdy - 2.0 * exp(x[2]) * area * lambdad3;

  /* upper triangle contribution */
  dvdx = (x[3] - x[2]) / hx;
  dvdy = (x[3] - x[1]) / hy;
  fquad += dvdx * dvdx + dvdy * dvdy;
  fexp += exp(x[1]) + exp(x[2]) + exp(x[3]);

  dvdx = 0.5 * dvdx * hy;
  dvdy = 0.5 * dvdy * hx;
  g[1] += -dvdy;
  g[2] += -dvdx;
  g[3] = dvdx + dvdy - exp(x[3]) * area * lambdad3;


  *f = area * (0.5 * fquad - lambdad3 * fexp);

  return 0;
} /* CombLocalFunctionGradient */


/*------- USER-DEFINED: routine to evaluate the Hessian
           at a local (rectangular element) level       -------*/

#undef __FUNCT__
#define __FUNCT__ "CombLocalHessian"
/*
  CombLocalHessian - Computes the Hessian of the local (partial) function
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
static int CombLocalHessian(PetscInt coor[2], double x[4], double H[4][4],void *ptr) {

  AppCtx *user = (AppCtx*)ptr;
  double hx, hy, lambdad3, area, dxdy, dydx;
  double diagxy, dexp, bandxy, bandyx;

  hx = user->hx;
  hy = user->hy;
  lambdad3 = user->param / 3.0;
  area = user->area;
  dxdy = hx/hy;
  dydx = hy/hx;
  diagxy = 0.5 * (dxdy + dydx);
  bandxy = -0.5 * dxdy;
  bandyx = -0.5 * dydx;

          /* Hessian contribution at 0,0 */
  dexp = exp(x[0]) * area * lambdad3;
  H[0][0] = diagxy - dexp;
  H[0][1] = H[1][0] = bandyx;
  H[0][2] = H[2][0] = bandxy;
  H[0][3] = H[3][0] = 0.0;

          /* Hessian contribution at 1,0 */
  dexp = exp(x[1]) * area * 2.0 * lambdad3;
  H[1][1] = diagxy - dexp;
  H[1][2] = H[2][1] = 0.0;
  H[1][3] = H[3][1] = bandxy;

          /* Hessian contribution at 0,1 */
  dexp = exp(x[2]) * area * 2.0 * lambdad3;
  H[2][2] = diagxy - dexp;
  H[2][3] = H[3][2] = bandyx;

          /* Hessian contribution at 1,1 */
  dexp = exp(x[3]) * area * lambdad3;
  H[3][3] = diagxy - dexp;

  return 0;
} /* CombLocalHessian */


/*------- USER-DEFINED: routine to evaluate the function 
          and gradient at the whole grid             -------*/

#undef __FUNCT__
#define __FUNCT__ "WholeCombFunctionGradient"
/*
  WholeCombFunctionGradient - Evaluates function and gradient over the 
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
static int WholeCombFunctionGradient(TAO_APPLICATION daapplication, DA da, Vec X, double *f, Vec G, void *ptr) {

  AppCtx *user = (AppCtx*)ptr;
  Vec localX, localG;
  int info;
  PetscInt  i, j;
  PetscInt xs, xm, gxs, gxm, xe, ys, ym, gys, gym, ye;
  double **x, **g;
  double floc = 0.0;
  PetscScalar zero = 0.0;

  double lambdad3, hx, hy, area;
  double fquad, fexp, dvdx, dvdy;

  lambdad3 = user->param / 3.0;
  hx = user->hx;
  hy = user->hy;
  area = user->area;

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
      fquad = dvdx * dvdx + dvdy * dvdy;
      fexp = exp(x[j][i]) + exp(x[j][i+1]) + exp(x[j+1][i]);

      dvdx = 0.5 * dvdx * hy;
      dvdy = 0.5 * dvdy * hx;
      g[j][i] += dvdx + dvdy - exp(x[j][i]) * area * lambdad3;
      g[j][i+1] += -dvdx - 2.0 * exp(x[j][i+1]) * area * lambdad3;
      g[j+1][i] += -dvdy - 2.0 * exp(x[j+1][i]) * area * lambdad3;

      /* upper triangle contribution */
      dvdx = (x[j+1][i+1] - x[j+1][i]) / hx;
      dvdy = (x[j+1][i+1] - x[j][i+1]) / hy;
      fquad += dvdx * dvdx + dvdy * dvdy;
      fexp += exp(x[j][i+1]) + exp(x[j+1][i]) + exp(x[j+1][i+1]);

      dvdx = 0.5 * dvdx * hy;
      dvdy = 0.5 * dvdy * hx;
      g[j][i+1] += -dvdy;
      g[j+1][i] += -dvdx;
      g[j+1][i+1] += dvdx + dvdy - exp(x[j+1][i+1]) * area * lambdad3;

      floc += area * (0.5 * fquad - fexp * lambdad3);

    }
  }

  info = MPI_Allreduce(&floc, f, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); CHKERRQ(info);

  info = DAVecRestoreArray(da, localX, (void**)&x); CHKERRQ(info);
  info = DAVecRestoreArray(da, localG, (void**)&g); CHKERRQ(info);

  info = DALocalToGlobalBegin(da, localG, G); CHKERRQ(info);
  info = DALocalToGlobalEnd(da, localG, G); CHKERRQ(info);

  info = DARestoreLocalVector(da, &localX); CHKERRQ(info);
  info = DARestoreLocalVector(da, &localG); CHKERRQ(info);

  info = PetscLogFlops((xe-xs) * (ye-ys) * 55 + 1); CHKERRQ(info);
  return 0;

} /* WholeCombFunctionGradient */


/*------- USER-DEFINED: routine to evaluate the Hessian 
          at the whole grid             -------*/
#undef __FUNCT__
#define __FUNCT__ "WholeCombHessian"
/*
  WholeCombHessian - Evaluates Hessian over the whole grid

  Input:
    daapplication - TAO application object
    da  - distributed array
    X   - the current point, at which the function and gradient are evaluated
    ptr - user-defined application context

  Output:
    H - Hessian at X
*/
static int WholeCombHessian(TAO_APPLICATION daapplication, DA da, Vec X, Mat H, void *ptr) {

  AppCtx *user = (AppCtx*)ptr;
  Vec localX;
  int info;
  PetscInt i, j, ind[4];
  PetscInt xs, xm, gxs, gxm, xe, ys, ym, gys, gym, ye;
  double smallH[4][4];
  double **x;

  double hx, hy, lambdad3, area, dxdy, dydx;
  double diagxy, dexp, bandxy, bandyx;
  PetscTruth assembled;


  hx = user->hx;
  hy = user->hy;
  lambdad3 = user->param / 3.0;
  area = user->area;
  dxdy = hx/hy;
  dydx = hy/hx;
  diagxy = 0.5 * (dxdy + dydx);
  bandxy = -0.5 * dxdy;
  bandyx = -0.5 * dydx;

  info = DAGetLocalVector(da, &localX); CHKERRQ(info);
  info = MatAssembled(H,&assembled); CHKERRQ(info);
  if (assembled){info = MatZeroEntries(H);  CHKERRQ(info);}

  info = DAGlobalToLocalBegin(da, X, INSERT_VALUES, localX); CHKERRQ(info);
  info = DAGlobalToLocalEnd(da, X, INSERT_VALUES, localX); CHKERRQ(info);

  info = DAVecGetArray(da, localX, (void**)&x); CHKERRQ(info);

  info = DAGetCorners(da, &xs, &ys, TAO_NULL, &xm, &ym, TAO_NULL); CHKERRQ(info);
  info = DAGetGhostCorners(da, &gxs, &gys, TAO_NULL, &gxm, &gym, TAO_NULL); CHKERRQ(info);

  xe = gxs + gxm - 1;
  ye = gys + gym - 1;
  for (j = ys; j < ye; j++) {
    for (i = xs; i < xe; i++) {

          /* Hessian contribution at 0,0 */
      dexp = exp(x[j][i]) * area * lambdad3;
      smallH[0][0] = diagxy - dexp;
      smallH[0][1] = smallH[1][0] = bandyx;
      smallH[0][2] = smallH[2][0] = bandxy;
      smallH[0][3] = smallH[3][0] = 0.0;

          /* Hessian contribution at 1,0 */
      dexp = exp(x[j][i+1]) * area * 2.0 * lambdad3;
      smallH[1][1] = diagxy - dexp;
      smallH[1][2] = smallH[2][1] = 0.0;
      smallH[1][3] = smallH[3][1] = bandxy;

          /* Hessian contribution at 0,1 */
      dexp = exp(x[j+1][i]) * area * 2.0 * lambdad3;
      smallH[2][2] = diagxy - dexp;
      smallH[2][3] = smallH[3][2] = bandyx;

          /* Hessian contribution at 1,1 */
      dexp = exp(x[j+1][i+1]) * area * lambdad3;
      smallH[3][3] = diagxy - dexp;

      ind[0] = (j-gys) * gxm + (i-gxs);
      ind[1] = ind[0] + 1;
      ind[2] = ind[0] + gxm;
      ind[3] = ind[2] + 1;
      info = MatSetValuesLocal(H,4,ind,4,ind,(PetscScalar*)smallH,ADD_VALUES); CHKERRQ(info);

    }
  }

  info = DAVecRestoreArray(da, localX, (void**)&x); CHKERRQ(info);

  info = MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  info = MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  info = MatSetOption(H, MAT_SYMMETRIC, PETSC_TRUE); CHKERRQ(info);

  info = DARestoreLocalVector(da, &localX); CHKERRQ(info);

  info = PetscLogFlops((xe-xs) * (ye-ys) * 14 + 7); CHKERRQ(info);
  return 0;

} /* WholeCombHessian */
