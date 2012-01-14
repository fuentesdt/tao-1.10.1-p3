/*$Id: eptorsion.c, v 1.1 2002/08/08 10:30 lopezca@mauddib.mcs.anl.gov $*/

/* Program usage: mpirun -np <proc> eptorsion [-help] [all TAO options] */

/*
  Include "tao.h" so we can use TAO solvers.
  petscda.h for distributed array
  ad_deriv.h for AD gradient
*/

#include "petscda.h"
#include "tao.h"
#include "taodaapplication.h"

static char help[] = "This example is based on the Elastic-Plastic Torsion (dept)\n\
problem from the MINPACK-2 test suite.\n\
The command line options are:\n\
  -mx <xg>, where <xg> = number of grid points in the 1st coordinate direction\n\
  -my <yg>, where <yg> = number of grid points in the 2nd coordinate direction\n\
  -nlevels <nlevels>, where <nlevels> = number of levels in multigrid\n\
  -byelement, if computation is made by functions on rectangular elements\n\
  -adic, if AD is used (AD is not used by default)\n\
  -u1 <u1>, where <u1> = upper limit in the 1st coordinate direction\n\
  -u2 <u2>, where <u2> = upper limit in the 2nd coordinate direction\n\
  -par <param>, where <param> = angle of twist per unit length\n\n";

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
   Routines: TaoSetOptions();
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
  PetscReal      param;          /* 'c' parameter */
  PetscReal      u1, u2;         /* domain upper limits (lower limits = 0) */
  double      hx, hy;        /* increment size in both directions */
  double      area;          /* area of the triangles */
  ADFGCtx     fgctx;         /* Used only when an ADIC generated gradient is used */
} AppCtx;
int ad_EPTorsLocalFunction(PetscInt[2], DERIV_TYPE[4], DERIV_TYPE*, void*);

/* User-defined routines found in this file */
static int AppCtxInitialize(void *ptr);
static int FormInitialGuess(DA, Vec);

static int EPTorsLocalFunctionGradient(PetscInt[2], double x[4], double *f, double g[4], void *ptr);
static int EPTorsLocalHessian(PetscInt[2], double x[4], double H[4][4], void *ptr);

static int WholeEPTorsFunctionGradient(TAO_APPLICATION,DA,Vec,double *,Vec,void*);
static int WholeEPTorsHessian(TAO_APPLICATION,DA,Vec,Mat,void*);

static int DASetBounds(TAO_APPLICATION, DA, Vec, Vec, void*);

static int MyGridMonitorBefore(TAO_APPLICATION, DA, PetscInt, void *);

#undef __FUNCT__
#define __FUNCT__ "main"

int main( int argc, char **argv ) {

  int             info;                           /* used to check for functions returning nonzeros */
  PetscInt             mx,my,Nx,Ny;
  double          ff,gnorm;
  PetscInt             iter, nlevels;                                                /* multigrid levels */
  DA              DAarray[20];
  Vec             X;
  PetscTruth      flg, PreLoad = PETSC_TRUE;                                               /* flags */
  TaoMethod       method = "tao_gpcg";                                       /* minimization method */
  AppCtx          user;                                                /* user-defined work context */
  TAO_SOLVER      tao;                                                 /* TAO_SOLVER solver context */
  TAO_APPLICATION EPTorsApp;                                            /* The PETSc application */
  TaoTerminateReason reason;
  KSP ksp;
  PC pc;

  /* Initialize TAO */
  PetscInitialize(&argc, &argv, (char *)0, help);
  TaoInitialize(&argc, &argv, (char *)0, help);

  PreLoadBegin(PreLoad,"Solve");
  
  info = AppCtxInitialize((void*)&user); CHKERRQ(info);

  nlevels=5;
  info = PetscOptionsGetInt(PETSC_NULL,"-nlevels",&nlevels,&flg); CHKERRQ(info);
  mx = my = 11;                               /* these correspond to 10 segments on each dimension */
  info = PetscOptionsGetInt(TAO_NULL, "-mx", &mx, &flg); CHKERRQ(info);
  info = PetscOptionsGetInt(TAO_NULL, "-my", &my, &flg); CHKERRQ(info);
  if (PreLoadIt == 0) {
    nlevels = 1; mx = 11; my = 11; }

  PetscPrintf(MPI_COMM_WORLD,"\n---- Elastic-Plastic Torsion Problem -----\n\n");

  /* Let PETSc determine the vector distribution */
  Nx = PETSC_DECIDE; Ny = PETSC_DECIDE;

  /* Create distributed array (DA) to manage parallel grid and vectors  */
  info = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,mx,
                    my,Nx,Ny,1,1,PETSC_NULL,PETSC_NULL,&DAarray[0]); CHKERRQ(info);
  for (iter=1;iter<nlevels;iter++){
    info = DARefine(DAarray[iter-1],PETSC_COMM_WORLD,&DAarray[iter]); CHKERRQ(info);
  }

  /* Create TAO solver and set desired solution method */
  info = TaoCreate(MPI_COMM_WORLD,method,&tao); CHKERRQ(info);
  info = TaoApplicationCreate(PETSC_COMM_WORLD,&EPTorsApp); CHKERRQ(info);
  info = TaoAppSetDAApp(EPTorsApp, DAarray, nlevels ); CHKERRQ(info);
  /* Sets routines for function, gradient and bounds evaluation */
  info = DAAppSetVariableBoundsRoutine(EPTorsApp,DASetBounds,(void *)&user); CHKERRQ(info);

  info = PetscOptionsHasName(TAO_NULL, "-byelement", &flg); CHKERRQ(info);
  if (flg) {

    /* Sets routines for function and gradient evaluation, element by element */
    info = PetscOptionsHasName(TAO_NULL, "-adic", &flg); CHKERRQ(info);
    if (flg) {
      info = DAAppSetADElementFunctionGradient(EPTorsApp,ad_EPTorsLocalFunction,192,(void *)&user.fgctx); CHKERRQ(info);
    } else {
      info = DAAppSetElementObjectiveAndGradientRoutine(EPTorsApp,EPTorsLocalFunctionGradient,42,(void *)&user); CHKERRQ(info);
    }
    /* Sets routines for Hessian evaluation, element by element */
    info = DAAppSetElementHessianRoutine(EPTorsApp,EPTorsLocalHessian,6,(void*)&user); CHKERRQ(info);

  } else {

    /* Sets routines for function and gradient evaluation, all in one routine */
    info = DAAppSetObjectiveAndGradientRoutine(EPTorsApp,WholeEPTorsFunctionGradient,(void *)&user); CHKERRQ(info);

    /* Sets routines for Hessian evaluation, all in one routine */
    info = DAAppSetHessianRoutine(EPTorsApp,WholeEPTorsHessian,(void*)&user); CHKERRQ(info);    
    
  }

  info = DAAppSetBeforeMonitor(EPTorsApp,MyGridMonitorBefore,(void*)&user); CHKERRQ(info);
  info = PetscOptionsHasName(TAO_NULL,"-tao_monitor", &flg); CHKERRQ(info);
  if (flg){
    info = DAAppPrintInterpolationError(EPTorsApp); CHKERRQ(info);
    info = DAAppPrintStageTimes(EPTorsApp); CHKERRQ(info);
  }
  info = TaoAppSetRelativeTolerance(EPTorsApp,1.0e-6); CHKERRQ(info);
  info = TaoSetTolerances(tao,0,0,0,0); CHKERRQ(info);
  info = TaoSetGradientTolerances(tao,0,0,0); CHKERRQ(info);

  info = TaoAppGetKSP(EPTorsApp,&ksp);CHKERRQ(info);
  info = KSPSetType(ksp,KSPCG); CHKERRQ(info);
  info = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,100);CHKERRQ(info);
  info = KSPGetPC(ksp,&pc);CHKERRQ(info);
  info = PCSetType(pc,PCBJACOBI);CHKERRQ(info);

  /* Check for any tao command line options */
  info = TaoSetOptions(EPTorsApp, tao); CHKERRQ(info);

  info = DAAppGetSolution(EPTorsApp,0,&X); CHKERRQ(info);
  info = FormInitialGuess(DAarray[0],X); CHKERRQ(info);
  info = DAAppSetInitialSolution(EPTorsApp,X); CHKERRQ(info);
  
  /* SOLVE THE APPLICATION */
  info = TaoDAAppSolve(EPTorsApp, tao);  CHKERRQ(info);

  /* Get information on termination */
  info = TaoGetSolutionStatus(tao,&iter,&ff,&gnorm,0,0,&reason); CHKERRQ(info);
  if (reason <= 0 ){
    PetscPrintf(MPI_COMM_WORLD,"Try a different TAO method, adjust some parameters, or check the function evaluation routines\n");
    PetscPrintf(MPI_COMM_WORLD," Iterations: %d,  Function Value: %4.2e, Residual: %4.2e \n",iter,ff,gnorm);
  }

  info = PetscOptionsHasName(PETSC_NULL,"-view_sol",&flg); CHKERRQ(info);
  if (flg){
    info = DAAppGetSolution(EPTorsApp,nlevels-1,&X); CHKERRQ(info);
    info=VecView(X,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(info);
  }

  /*  To View TAO solver information */
  // info = TaoView(tao); CHKERRQ(info);

  /* Free TAO data structures */
  info = TaoDestroy(tao); CHKERRQ(info);
  info = TaoAppDestroy(EPTorsApp); CHKERRQ(info);

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
  help diplay info of iterations at every grid level
*/
#undef __FUNCT__
#define __FUNCT__ "MyGridMonitorBefore"
static int MyGridMonitorBefore(TAO_APPLICATION myapp, DA da, PetscInt level, void *ctx) {

  AppCtx *user = (AppCtx*)ctx;
  PetscInt mx,my;
  int info;

  info = DAGetInfo(da,PETSC_NULL,&mx,&my,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,
		   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(info);
  user->hx = user->u1 / (mx - 1);
  user->hy = user->u2 / (my - 1);
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
  PetscTruth    flg;            /* flag for PETSc calls */
  int info;

  /* Specify default parameters */
  user->param = 25.0;
  user->u1 = user->u2 = 1.0;

  /* Check for command line arguments that override defaults */
  info = PetscOptionsGetReal(TAO_NULL, "-par", &user->param, &flg); CHKERRQ(info);
  info = PetscOptionsGetReal(TAO_NULL, "-u1", &user->u1, &flg); CHKERRQ(info);
  info = PetscOptionsGetReal(TAO_NULL, "-u2", &user->u2, &flg); CHKERRQ(info);

  return 0;
} /* AppCtxInitialize */

#undef __FUNCT__
#define __FUNCT__ "FormInitialGuess"
static int FormInitialGuess(DA da, Vec X)
{
  int info;
  PetscInt    i, j, mx, my;
  PetscInt    xs, ys, xm, ym, xe, ye;
  PetscReal hx, hy, temp, val;
  double **x;

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
      val = PetscMin((PetscMin(i+1,mx-i))*hx,temp);
      x[j][i] = val;
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
static int DASetBounds(TAO_APPLICATION daapplication, DA da, Vec XL, Vec XU, void *ptr)
{
  AppCtx *user = (AppCtx*)ptr;
  int info;
  PetscInt i, j, xs, xm, ys, ym;
  double hx, hy, u1, u2, dist, d1, d2, hd, vd;
  double **xl, **xu;

  hx = user->hx;
  hy = user->hy;
  u1 = user->u1;
  u2 = user->u2;

  info = DAVecGetArray(da, XL, (void**)&xl); CHKERRQ(info);
  info = DAVecGetArray(da, XU, (void**)&xu); CHKERRQ(info);
  info = DAGetCorners(da, &xs, &ys, TAO_NULL, &xm, &ym, TAO_NULL); CHKERRQ(info);

  for (j = ys; j < ys+ym; j++){
    for (i = xs; i < xs+xm; i++){
      d1 = i * hx; d2 = u1 - d1; hd = PetscMin(d1,d2);
      d1 = j * hy; d2 = u2 - d1; vd = PetscMin(d1,d2);
      dist = PetscMin(hd,vd);
      xl[j][i] = -dist;
      xu[j][i] = dist;
    }
  }

  info = DAVecRestoreArray(da, XL, (void**)&xl); CHKERRQ(info);
  info = DAVecRestoreArray(da, XU, (void**)&xu); CHKERRQ(info);

  info = PetscLogFlops(xm * ym * 4); CHKERRQ(info);
  return 0;

} /* DASetBounds */


#undef __FUNCT__
#define __FUNCT__ "EPTorsLocalFunctionGradient"
/*
  EPTorsLocalFunctionGradient - Evaluates function and gradient over the 
      local rectangular element

  Input:
    coor - vector with the indices of the position of current element
             in the first, second and third directions
    x - current point (values over the current rectangular element)
    df - degrees of freedom at each point
    ptr - user-defined application context

  Output:
    f - value of the objective funtion at the local rectangular element
    g - gradient of the local function
*/
static int EPTorsLocalFunctionGradient(PetscInt coor[2], double x[4], double *f, double g[4], void *ptr) {

  AppCtx *user = (AppCtx*)ptr;

  double fquad, flin;
  double hx, hy, dvdx, dvdy, area;
  double cdiv3, cnt;

  cdiv3 = user->param / 3.0;
  hx = user->hx;
  hy = user->hy;
  area = user->area;
  cnt = area * cdiv3;

  /* lower triangle contribution */
  dvdx = (x[0] - x[1]) / hx;
  dvdy = (x[0] - x[2]) / hy;
  fquad = dvdx * dvdx + dvdy * dvdy;
  flin = x[0] + x[1] + x[2];

  dvdx = 0.5 * dvdx * hy;
  dvdy = 0.5 * dvdy * hx;
  g[0] = dvdx + dvdy - cnt;
  g[1] = -dvdx - 2.0 * cnt;
  g[2] = -dvdy - 2.0 * cnt;

  /* upper triangle contribution */
  dvdx = (x[3] - x[2]) / hx;
  dvdy = (x[3] - x[1]) / hy;
  fquad += dvdx * dvdx + dvdy * dvdy;
  flin += x[1] + x[2] + x[3];

  dvdx = 0.5 * dvdx * hy;
  dvdy = 0.5 * dvdy * hx;
  g[1] += -dvdy;
  g[2] += -dvdx;
  g[3] = dvdx + dvdy - cnt;

  *f = area * (0.5 * fquad - flin * cdiv3);

  return 0;
} /* EPTorsLocalFunctionGradient */



/*------- USER-DEFINED: routine to evaluate the Hessian
           at a local (rectangular element) level       -------*/
#undef __FUNCT__
#define __FUNCT__ "EPTorsLocalHessian"
/*
  EPTorsLocalHessian - Computes the Hessian of the local (partial) function
         defined over the current rectangle

  Input:
    coor - vector with the indices of the position of current element
             in the first, second and third directions
    x - current local solution (over the rectangle only)
    df - degrees of freedom at each point
    ptr - user-defined application context

  Output:
    H - Hessian matrix of the local function (wrt the four
           points of the rectangle only)
*/
static int EPTorsLocalHessian(PetscInt coor[2], double x[4], double H[4][4], void *ptr) {

  AppCtx *user = (AppCtx*)ptr;
  double hx, hy, dxdy, dydx;
  double diagxy, bandxy, bandyx;

  hx = user->hx;
  hy = user->hy;
  dxdy = hx/hy;
  dydx = hy/hx;
  diagxy = 0.5 * (dxdy + dydx);
  bandxy = -0.5 * dxdy;
  bandyx = -0.5 * dydx;

          /* Hessian contribution at 0,0 */
  H[0][0] = diagxy;
  H[0][1] =  H[1][0] = bandyx;
  H[0][2] =  H[2][0] = bandxy;
  H[0][3] =  H[3][0] = 0.0;

          /* Hessian contribution at 1,0 */
  H[1][1] = diagxy;
  H[1][2] =  H[2][1] = 0.0;
  H[1][3] =  H[3][1] = bandxy;

          /* Hessian contribution at 0,1 */
  H[2][2] = diagxy;
  H[2][3] =  H[3][2] = bandyx;

          /* Hessian contribution at 1,1 */
  H[3][3] = diagxy;

  return 0;

} /* EPTorsLocalHessian */


/*------- USER-DEFINED: routine to evaluate the function 
          and gradient at the whole grid             -------*/
#undef __FUNCT__
#define __FUNCT__ "WholeEPTorsFunctionGradient"
/*
  WholeEPTorsFunctionGradient - Evaluates function and gradient over the 
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
static int WholeEPTorsFunctionGradient(TAO_APPLICATION daapplication, DA da, Vec X, double *f, Vec G, void *ptr) {

  AppCtx *user = (AppCtx*)ptr;
  Vec localX, localG;
  int info;
  PetscInt i, j;
  PetscInt xs, xm, gxs, gxm, xe, ys, ym, gys, gym, ye;
  double **x, **g;
  double floc = 0.0;
  PetscScalar zero = 0.0;

  double fquad, flin;
  double hx, hy, dvdx, dvdy, area;
  double cdiv3, cnt;

  cdiv3 = user->param / 3.0;
  hx = user->hx;
  hy = user->hy;
  area = user->area;
  cnt = area * cdiv3;

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
      flin = x[j][i] + x[j][i+1] + x[j+1][i];

      dvdx = 0.5 * dvdx * hy;
      dvdy = 0.5 * dvdy * hx;
      g[j][i] += dvdx + dvdy - cnt;
      g[j][i+1] += -dvdx - 2.0 * cnt;
      g[j+1][i] += -dvdy - 2.0 * cnt;

      /* upper triangle contribution */
      dvdx = (x[j+1][i+1] - x[j+1][i]) / hx;
      dvdy = (x[j+1][i+1] - x[j][i+1]) / hy;
      fquad += dvdx * dvdx + dvdy * dvdy;
      flin += x[j][i+1] + x[j+1][i] + x[j+1][i+1];

      dvdx = 0.5 * dvdx * hy;
      dvdy = 0.5 * dvdy * hx;
      g[j][i+1] += -dvdy;
      g[j+1][i] += -dvdx;
      g[j+1][i+1] += dvdx + dvdy - cnt;

      floc += area * (0.5 * fquad - flin * cdiv3);

    }
  }

  info = MPI_Allreduce(&floc, f, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); CHKERRQ(info);

  info = DAVecRestoreArray(da, localX, (void**)&x); CHKERRQ(info);
  info = DAVecRestoreArray(da, localG, (void**)&g); CHKERRQ(info);

  info = DALocalToGlobalBegin(da, localG, G); CHKERRQ(info);
  info = DALocalToGlobalEnd(da, localG, G); CHKERRQ(info);

  info = DARestoreLocalVector(da, &localX); CHKERRQ(info);
  info = DARestoreLocalVector(da, &localG); CHKERRQ(info);

  info = PetscLogFlops((xe-xs) * (ye-ys) * 47 + 2); CHKERRQ(info);
  return 0;
} /* WholeEPTorsFunctionGradient */


/*------- USER-DEFINED: routine to evaluate the Hessian 
          at the whole grid             -------*/
#undef __FUNCT__
#define __FUNCT__ "WholeEPTorsHessian"
/*
  WholeEPTorsHessian - Evaluates Hessian over the whole grid

  Input:
    daapplication - TAO application object
    da  - distributed array
    X   - the current point, at which the function and gradient are evaluated
    ptr - user-defined application context

  Output:
    H - Hessian at X
*/
static int WholeEPTorsHessian(TAO_APPLICATION daapplication, DA da, Vec X, Mat H, void *ptr) {

  AppCtx *user = (AppCtx*)ptr;
  int info;
  PetscInt i, j, ind[4];
  PetscInt xs, xm, gxs, gxm, xe, ys, ym, gys, gym, ye;
  double smallH[4][4];

  double hx, hy, dxdy, dydx;
  double diagxy, bandxy, bandyx;
  PetscTruth assembled;

  hx = user->hx;
  hy = user->hy;
  dxdy = hx/hy;
  dydx = hy/hx;
  diagxy = 0.5 * (dxdy + dydx);
  bandxy = -0.5 * dxdy;
  bandyx = -0.5 * dydx;

  info = MatAssembled(H,&assembled); CHKERRQ(info);
  if (assembled){info = MatZeroEntries(H);  CHKERRQ(info);}


  info = DAGetCorners(da, &xs, &ys, TAO_NULL, &xm, &ym, TAO_NULL); CHKERRQ(info);
  info = DAGetGhostCorners(da, &gxs, &gys, TAO_NULL, &gxm, &gym, TAO_NULL); CHKERRQ(info);

  xe = gxs + gxm - 1;
  ye = gys + gym - 1;
  for (j = ys; j < ye; j++) {
    for (i = xs; i < xe; i++) {

          /* Hessian contribution at 0,0 */
      smallH[0][0] = diagxy;
      smallH[0][1] = smallH[1][0] = bandyx;
      smallH[0][2] = smallH[2][0] = bandxy;
      smallH[0][3] = smallH[3][0] = 0.0;

          /* Hessian contribution at 1,0 */
      smallH[1][1] = diagxy;
      smallH[1][2] = smallH[2][1] = 0.0;
      smallH[1][3] = smallH[3][1] = bandxy;

          /* Hessian contribution at 0,1 */
      smallH[2][2] = diagxy;
      smallH[2][3] = smallH[3][2] = bandyx;

          /* Hessian contribution at 1,1 */
      smallH[3][3] = diagxy;

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


  info = PetscLogFlops(6); CHKERRQ(info);
  return 0;

} /* WholeEPTorsHessian */
