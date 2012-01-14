/*$Id: minsurf3.c,v 1.1 2002/05/13 19:32:57 benson Exp $*/

/* Program usage: mpirun -np <proc> minsurf3 [-help] [all TAO options] */

/* minsurf1: boundary values like in COPS example:
      - parabola for top and bottom boundaries
      - zero o.w.  */

/*
  Include "tao.h" so we can use TAO solvers.
  petscda.h for distributed array
  ad_deriv.h for AD gradient
*/

#include "petscda.h"
#include "tao.h"
#include "taodaapplication.h"

static char  help[] = 
"Given a rectangular 2-D domain and boundary values along \n\
the edges of the domain, the objective is to find the surface with \n\
the minimal area that satisfies the boundary conditions.\n\
The command line options are:\n\
  -mx <xg>, where <xg> = number of grid points in the 1st coordinate direction\n\
  -my <yg>, where <yg> = number of grid points in the 2nd coordinate direction\n\
  -nlevels <nlevels>, where <nlevels> = number of levels in multigrid\n\
  -byelement, if computation is made by functions on rectangular elements\n\
  -adic, if AD is used (AD is not used by default)\n\
  -bottom <b>, where <b> = bottom bound on the rectangular domain\n\
  -top <t>, where <t> = bottom bound on the rectangular domain\n\
  -left <l>, where <l> = left bound on the rectangular domain\n\
  -right <r>, where <r> = right bound on the rectangular domain\n\
  -cops <r>, if COPS boundaries should be use (default MINPACK bounds)\n\
  -obst <obst>, where <obst> = 1 is obstacle is present, zero otherwise \n\n";

/*T
   Concepts: TAO - Solving a bounded minimization problem
   Routines: TaoInitialize(); TaoFinalize();
   Routines: TaoCreate(); TaoDestroy();
   Routines: DAAppSetVariableBoundsRoutine();
   Routines: DAAppSetElementObjectiveAndGradientRoutine();
   Routines: DAAppSetElementHessianRoutine();
   Routines: DAAppSetObjectiveAndGradientRoutine();
   Routines: DAAppSetADElementFunctionGradient();
   Routines: DAAppSetHessianRoutine();
   Routines: TaoSetOptions();
   Routines: TaoAppGetSolutionStatus(); TaoDAAppSolve();
   Routines: DAAppSetBeforeMonitor(); TaoView();
   Routines: DAAppGetSolution();
   Routines: DAAppGetInterpolationMatrix();
   Processors: n
T*/

/*
   User-defined application context - contains data needed by the
   application-provided call-back routines.
*/


/*  
    This structure is used only when an ADIC generated gradient is used.
    An InactiveDouble type is a double 
*/
typedef struct {
  
  InactiveDouble      hx, hy;        /* increment size in both directions */
  InactiveDouble      area;          /* area of the triangles */

} ADFGCtx;
int ad_MSurfLocalFunction(PetscInt[2], DERIV_TYPE[4], DERIV_TYPE*, void*);

typedef struct {
  PetscReal      b, t, l, r;    /* domain boundaries */
  PetscReal      bheight;       /* height of obstacle    */
  PetscReal      fx, fy;        /* relative size of obstacle */
  double      hx, hy;        /* increment size in both directions */
  double      area;          /* area of the triangles */
  PetscInt         mx, my;        /* discretization including boundaries */
  PetscInt         bmx, bmy;      /* discretization of obstacle */
  ADFGCtx     fgctx;         /* Used only when an ADIC generated gradient is used */
} AppCtx;

/* User-defined routines */
static int AppCtxInitialize(void *);

static int MSurfLocalFunctionGradient(PetscInt[2], double[4], double *, double[4], void *);
static int WholeMSurfFunctionGradient(TAO_APPLICATION,DA,Vec,double *,Vec,void*);

static int MSurfLocalHessian(PetscInt[2], double[4], double[4][4], void *);
static int WholeMSurfHessian(TAO_APPLICATION,DA,Vec,Mat,void*);

static int COPS_Bounds(TAO_APPLICATION, DA, Vec, Vec, void*);
static int MINPACK_Bounds(TAO_APPLICATION, DA, Vec, Vec, void *);

static int MyGridMonitorBefore(TAO_APPLICATION, DA, PetscInt, void *);

#undef __FUNCT__
#define __FUNCT__ "main"

int main( int argc, char **argv ) {

  int                  info;         /* used to check for functions returning nonzeros */
  PetscInt                  iter,nlevels;
  PetscInt                  Nx,Ny;
  double               ff,gnorm;
  DA                   DAarray[20];
  Vec                  X;
  PetscTruth           flg, PreLoad = PETSC_TRUE;                             /* flags */
  TaoMethod            method = "tao_bnls";                     /* minimization method */
  TaoTerminateReason   reason;
  TAO_SOLVER           tao;                               /* TAO_SOLVER solver context */
  TAO_APPLICATION   MSurfApp;                              /* The PETSc application */
  AppCtx               user;                              /* user-defined work context */
  KSP ksp;
  PC  pc;

  /* Initialize TAO */
  PetscInitialize(&argc, &argv, (char *)0, help);
  TaoInitialize(&argc, &argv, (char *)0, help);

  PreLoadBegin(PreLoad,"Solve");

  info = AppCtxInitialize((void*)&user); CHKERRQ(info);

  nlevels=4;
  info = PetscOptionsGetInt(PETSC_NULL,"-nlevels",&nlevels,&flg); CHKERRQ(info);
  if (PreLoadIt == 0) {
    nlevels = 1; user.mx = 11; user.my = 11;
  }

  PetscPrintf(MPI_COMM_WORLD,"\n---- Minimal Surface Problem (simple boundary) -----\n\n");

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
  info = TaoApplicationCreate(PETSC_COMM_WORLD, &MSurfApp); CHKERRQ(info);
  info = TaoAppSetDAApp(MSurfApp,DAarray,nlevels); CHKERRQ(info);

  /* Sets routines for function, gradient and bounds evaluation */
  info = PetscOptionsHasName(TAO_NULL, "-cops", &flg); CHKERRQ(info);
  if (flg){
    info = DAAppSetVariableBoundsRoutine(MSurfApp,COPS_Bounds,(void *)&user); CHKERRQ(info);
  } else {
    info = DAAppSetVariableBoundsRoutine(MSurfApp,MINPACK_Bounds,(void *)&user); CHKERRQ(info);
  }

  info = PetscOptionsHasName(TAO_NULL, "-byelement", &flg); CHKERRQ(info);
  if (flg) {

    /* Sets routines for function and gradient evaluation, element by element */
    info = PetscOptionsHasName(TAO_NULL, "-adic", &flg); CHKERRQ(info);
    if (flg) {
      info = DAAppSetADElementFunctionGradient(MSurfApp,ad_MSurfLocalFunction,150,(void *)&user.fgctx); CHKERRQ(info);
    } else {
      info = DAAppSetElementObjectiveAndGradientRoutine(MSurfApp,MSurfLocalFunctionGradient,36,(void *)&user); CHKERRQ(info);
    }

    /* Sets routines for Hessian evaluation, element by element */
    info = DAAppSetElementHessianRoutine(MSurfApp,MSurfLocalHessian,87,(void*)&user); CHKERRQ(info);

  } else {

    /* Sets routines for function and gradient evaluation, all in one routine */
    info = DAAppSetObjectiveAndGradientRoutine(MSurfApp,WholeMSurfFunctionGradient,(void *)&user); CHKERRQ(info);

    /* Sets routines for Hessian evaluation, all in one routine */
    info = DAAppSetHessianRoutine(MSurfApp,WholeMSurfHessian,(void*)&user); CHKERRQ(info);    
    
  }

  info = DAAppSetBeforeMonitor(MSurfApp,MyGridMonitorBefore,(void*)&user); CHKERRQ(info);
  info = PetscOptionsHasName(TAO_NULL,"-tao_monitor", &flg); CHKERRQ(info);
  if (flg){
    info = DAAppPrintInterpolationError(MSurfApp); CHKERRQ(info);
    info = DAAppPrintStageTimes(MSurfApp); CHKERRQ(info);
  }

  info = TaoAppSetRelativeTolerance(MSurfApp,1.0e-6); CHKERRQ(info);
  info = TaoSetTolerances(tao,0,0,0,0); CHKERRQ(info);
  info = TaoSetGradientTolerances(tao,0,0,0); CHKERRQ(info);

  info = TaoAppGetKSP(MSurfApp,&ksp);CHKERRQ(info);
  info = KSPSetType(ksp,KSPCG); CHKERRQ(info);
  info = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,100);CHKERRQ(info);
  info = KSPGetPC(ksp,&pc);CHKERRQ(info);
  info = PCSetType(pc,PCBJACOBI);CHKERRQ(info);

  /* Check for any tao command line options */
  info = TaoSetOptions(MSurfApp,tao); CHKERRQ(info);

  /* SOLVE THE APPLICATION */
  info = TaoDAAppSolve(MSurfApp,tao);  CHKERRQ(info);

  /* Get information on termination */
  info = TaoGetSolutionStatus(tao,&iter,&ff,&gnorm,0,0,&reason); CHKERRQ(info);
  if (reason <= 0 ){
    PetscPrintf(MPI_COMM_WORLD,"Try a different TAO method, adjust some parameters, or check the function evaluation routines\n");
    PetscPrintf(MPI_COMM_WORLD," Iterations: %d,  Function Value: %4.2e, Residual: %4.2e \n",iter,ff,gnorm);
  }

  info = PetscOptionsHasName(PETSC_NULL,"-view_sol",&flg); CHKERRQ(info);
  if (flg){
    info = DAAppGetSolution(MSurfApp,nlevels-1,&X); CHKERRQ(info);
    info=VecView(X,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(info);
  }

  /*  To View TAO solver information */
  // info = TaoView(tao); CHKERRQ(info);

  /* Free TAO data structures */
  info = TaoDestroy(tao); CHKERRQ(info);
  info = TaoAppDestroy(MSurfApp); CHKERRQ(info);

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
  int info;
  PetscInt mx,my;

  info = DAGetInfo(da,PETSC_NULL,&mx,&my,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,
		   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(info);
  user->mx = mx;
  user->my = my;
  user->hx = (user->r - user->l) / (user->mx - 1);
  user->hy = (user->t - user->b) / (user->my - 1);
  user->area = 0.5 * user->hx * user->hy;
  user->fgctx.hx   = user->hx;
  user->fgctx.hy   = user->hy;
  user->fgctx.area = user->area;

  user->bmx=(PetscInt)((mx+1)*user->fx); user->bmy=(PetscInt)((my+1)*user->fy); 

  PetscPrintf(MPI_COMM_WORLD,"Grid: %d,    mx: %d     my: %d   \n",level,mx,my);

  return 0;
}

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
  PetscTruth flg;
  int info;

  /* Specify default parameters */
  user->mx = user->my = 11;
  user->b = -0.5;
  user->t = 0.5;
  user->l = -0.5;
  user->r = 0.5;
  user->fx=0.5;
  user->fy=0.5;
  user->bheight=0.0;

  /* Check for command line arguments that override defaults */
  info = PetscOptionsGetInt(TAO_NULL, "-mx", &user->mx, &flg); CHKERRQ(info);
  info = PetscOptionsGetInt(TAO_NULL, "-my", &user->my, &flg); CHKERRQ(info);
  info = PetscOptionsGetReal(TAO_NULL, "-bottom", &user->b, &flg); CHKERRQ(info);
  info = PetscOptionsGetReal(TAO_NULL, "-top", &user->t, &flg); CHKERRQ(info);
  info = PetscOptionsGetReal(TAO_NULL, "-left", &user->l, &flg); CHKERRQ(info);
  info = PetscOptionsGetReal(TAO_NULL, "-right", &user->r, &flg); CHKERRQ(info);
  info = PetscOptionsGetReal(PETSC_NULL,"-bmx",&user->fx,&flg); CHKERRQ(info);
  info = PetscOptionsGetReal(PETSC_NULL,"-bmy",&user->fy,&flg); CHKERRQ(info);
  info = PetscOptionsGetReal(PETSC_NULL,"-bheight",&user->bheight,&flg); CHKERRQ(info);

  user->hx = (user->r - user->l) / (user->mx - 1);
  user->hy = (user->t - user->b) / (user->my - 1);
  user->area = 0.5 * user->hx * user->hy;
  info = PetscLogFlops(8); CHKERRQ(info);

  return 0;

} /* AppCtxInitialize */


/*------- USER-DEFINED: set the upper and lower bounds for the variables  -------*/
#undef __FUNCT__
#define __FUNCT__ "COPS_Bounds"

/*
  FormBounds - Forms bounds on the variables

  Input:
    user - user-defined application context

  Output:
    XL - vector of lower bounds
    XU - vector of upper bounds

*/
static int COPS_Bounds(TAO_APPLICATION tao, DA da, Vec XL, Vec XU, void *ptr) {

  AppCtx *user = (AppCtx*)ptr;
  PetscTruth flg;
  int info;
  PetscInt i, j, mx, my, xs, xm, ys, ym;
  double lb = -TAO_INFINITY;
  double ub = TAO_INFINITY;
  double **xl, **xu;
  double xi, xi1, xi2;
  double cx, cy, radx, rady, height = 1.0;

  info = DAGetInfo(da,PETSC_NULL,&mx,&my,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,
		   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(info);
  user->hx = (user->r - user->l) / (mx - 1);
  user->hy = (user->t - user->b) / (my - 1);
  user->area = 0.5 * user->hx * user->hy;

  info = DAVecGetArray(da, XL, (void**)&xl); CHKERRQ(info);
  info = DAVecGetArray(da, XU, (void**)&xu); CHKERRQ(info);
  info = DAGetCorners(da, &xs, &ys, TAO_NULL, &xm, &ym, TAO_NULL); CHKERRQ(info);

  /* Compute default bounds */
  for (j = ys; j < ys+ym; j++){
    for (i = xs; i < xs+xm; i++){

      if (j == 0 || j == my - 1) {
        xi = user->l + i * user->hx;
        xl[j][i] = xu[j][i] =  -4 * (xi - user->l) * (xi - user->r);
      } else if (i == 0 || i == mx - 1) {
        xl[j][i] = xu[j][i] = 0.0;
      } else {
        xl[j][i] = lb;
        xu[j][i] = ub;
      }

    }
  }

  /* Adjust lower bound if obstacle is present */
  info = PetscOptionsHasName(PETSC_NULL, "-obst", &flg); CHKERRQ(info);
  if (flg) {
    radx = (user->r - user->l) * 0.25;
    cx = user->l + 2.0 * radx;
    rady = (user->t - user->b) * 0.25;
    cy = user->b + 2.0 * rady;
    for (j = ys; j < ys+ym; j++){
      for (i = xs; i < xs+xm; i++){
        xi1 = user->l + i * user->hx;
        xi2 = user->b + j * user->hy;
        if ( fabs(xi1 - cx) <= radx && fabs(xi2 - cy) <= rady ) {
          xl[j][i] = height;
        }
      }
    }
    info = PetscLogFlops(8 + xm * ym * 6); CHKERRQ(info);
  }

  info = DAVecRestoreArray(da, XL, (void**)&xl); CHKERRQ(info);
  info = DAVecRestoreArray(da, XU, (void**)&xu); CHKERRQ(info);

  info = PetscLogFlops(12 * ym + 6); CHKERRQ(info);

  return 0;

} /* DAGetBounds2d */


#undef __FUNCT__
#define __FUNCT__ "MINPACK_Bounds"
static int MINPACK_Bounds(TAO_APPLICATION tao, DA da, Vec XL,Vec XU, void *user1){

  AppCtx *user=(AppCtx*)user1;
  int info;
  PetscInt i,j,k,limit=0,maxits=5;
  PetscInt xs,ys,xm,ym;
  PetscInt mx, my, bmy, bmx;
  PetscInt row=0, bsize=0, lsize=0, tsize=0, rsize=0;
  double bheight=user->bheight;
  double one=1.0, two=2.0, three=3.0, tol=1e-10;
  double fnorm,det,hx,hy,xt=0,yt=0;
  double u1,u2,nf1,nf2,njac11,njac12,njac21,njac22;
  double b=-0.5, t=0.5, l=-0.5, r=0.5;
  PetscScalar scl,*xl,*xu, lb=TAO_NINFINITY, ub=TAO_INFINITY;
  PetscTruth flg;
  double **xll;

  info = DAGetInfo(da,PETSC_NULL,&mx,&my,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,
		   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(info);

  bheight=user->bheight,lb=-1000; ub=1000;
  bmx=user->bmx; bmy=user->bmy;

  info = DAGetCorners(da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL); CHKERRQ(info);

  info = VecSet(XL, lb); CHKERRQ(info);
  info = VecSet(XU, ub); CHKERRQ(info);

  /*
    Get pointers to vector data
  */
  info = DAVecGetArray(da,XL,(void**)&xll); CHKERRQ(info);

  if (ys==0) bsize=xm;
  if (xs==0) lsize=ym;
  if (xs+xm==mx) rsize=ym;
  if (ys+ym==my) tsize=xm;
  
  hx= (r-l)/(mx-1); hy=(t-b)/(my-1);

  /* Compute the optional lower box */
  for (i=xs; i< xs+xm; i++){    
    for (j=ys; j<ys+ym; j++){
      
      if (i>=(mx-bmx)/2 && i<mx-(mx-bmx)/2 && j>=(my-bmy)/2 && j<my-(my-bmy)/2 ){
        xll[j][i] = bheight;
      }

    }
  }

  info = DAVecRestoreArray(da,XL,(void**)&xll); CHKERRQ(info);

  /* Boundary Values */
  info = VecGetArray(XL,&xl); CHKERRQ(info);
  info = VecGetArray(XU,&xu); CHKERRQ(info);

  for (j=0; j<4; j++){
    if (j==0){
      yt=b;
      xt=l+hx*xs;
      limit=bsize;
      scl=1.0;
      info = PetscOptionsGetReal(PETSC_NULL,"-bottom",&scl,&flg); 
    } else if (j==1){
      yt=t;
      xt=l+hx*xs;
      limit=tsize;
      scl=1.0;
      info = PetscOptionsGetReal(PETSC_NULL,"-top",&scl,&flg); 
    } else if (j==2){
      yt=b+hy*ys;
      xt=l;
      limit=lsize;
      scl=1.0;
      info = PetscOptionsGetReal(PETSC_NULL,"-left",&scl,&flg); 
    } else if (j==3){
      yt=b+hy*ys;
      xt=r;
      limit=rsize;
      scl=1.0;
      info = PetscOptionsGetReal(PETSC_NULL,"-right",&scl,&flg); 
    }

    for (i=0; i<limit; i++){
      u1=xt;
      u2=-yt;
      for (k=0; k<maxits; k++){
	nf1=u1 + u1*u2*u2 - u1*u1*u1/three-xt;
	nf2=-u2 - u1*u1*u2 + u2*u2*u2/three-yt;
	fnorm=sqrt(nf1*nf1+nf2*nf2);
	if (fnorm <= tol) break;
	njac11=one+u2*u2-u1*u1;
	njac12=two*u1*u2;
	njac21=-two*u1*u2;
	njac22=-one - u1*u1 + u2*u2;
	det = njac11*njac22-njac21*njac12;
	u1 = u1-(njac22*nf1-njac12*nf2)/det;
	u2 = u2-(njac11*nf2-njac21*nf1)/det;
      }

      if (j==0){
	row = i;
      } else if (j==1){
	row = (ym-1)*xm+i;
      } else if (j==2){
	row = (xm*i);
      } else if (j==3){
	row = xm*(i+1)-1;
      }
      
      xl[row]=(u1*u1-u2*u2)*scl;
      xu[row]=(u1*u1-u2*u2)*scl;

      if (j==0 || j==1) {
	xt=xt+hx;
      } else if (j==2 || j==3){
	yt=yt+hy;
      }
      
    }
    
  }

  info = VecRestoreArray(XL,&xl); CHKERRQ(info);
  info = VecRestoreArray(XU,&xu); CHKERRQ(info);

  return 0;
}

/*------- USER-DEFINED: routine to evaluate the function and gradient
           at a local (rectangular element) level              -------*/
#undef __FUNCT__
#define __FUNCT__ "MSurfLocalFunctionGradient"

/*
  MSurfLocalFunctionGradient - Evaluates function and gradient over the 
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
static int MSurfLocalFunctionGradient(PetscInt coor[2], double x[4], double *f, double g[4], void *ptr) {

  AppCtx *user = (AppCtx*)ptr;

  double hx, hy, area;
  double dvdx, dvdy, flow, fup;
  double d1,d2;

  hx = user->hx;
  hy = user->hy;
  area = user->area;

  /* lower triangle contribution */
  dvdx = (x[0] - x[1]);
  dvdy = (x[0] - x[2]);
  g[1] = (dvdx * (hy/hx))/(-2);
  g[2] = (dvdy * (hx/hy))/(-2);
  dvdx /= hx;
  dvdy /= hy;
  flow = sqrt( 1 + dvdx * dvdx + dvdy * dvdy );
  *f=flow;
  g[1] /= flow;
  g[2] /= flow;
  g[0] = -(g[1]+g[2]);

  /* upper triangle contribution */
  dvdx = (x[3] - x[2]);
  dvdy = (x[3] - x[1]);
  d1 = (dvdy*(hx/hy))/(-2);
  d2 = (dvdx*(hy/hx))/(-2);
  dvdx /= hx;
  dvdy /= hy;
  fup = sqrt( 1 + dvdx * dvdx + dvdy * dvdy );
  *f += fup;
  g[1] += d1/fup;
  g[2] += d2/fup; 
  g[3] = -(d1+d2)/fup;

  *f *= area;

  return 0;
} /* MSurfLocalFunctionGradient */



#undef __FUNCT__
#define __FUNCT__ "MSurfLocalHessian"
/*
  MSurfLocalHessian - Computes the Hessian of the local (partial) function
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
static int MSurfLocalHessian(PetscInt coor[2], double x[4], double H[4][4],void *ptr) {

  AppCtx *user = (AppCtx*)ptr;
  double hx, hy, area, byhxhx, byhyhy;
  double dvdx, dvdy, flow, fup;
  double areadivf, areadivf3;

  hx = user->hx;
  hy = user->hy;
  area = user->area;
  
  byhxhx = 1.0 / (hx * hx);
  byhyhy = 1.0 / (hy * hy);

  /* 0 is 0,0; 1 is 1,0; 2 is 0,1; 3 is 1,1 */
  dvdx = (x[0] - x[1]) / hx;  /* lower triangle contribution */
  dvdy = (x[0] - x[2]) / hy;
  flow = sqrt( 1 + dvdx * dvdx + dvdy * dvdy );
  dvdx = dvdx / hx;
  dvdy = dvdy / hy;
  areadivf = area / flow;
  areadivf3 = areadivf / (flow * flow);
  H[0][0] = areadivf * (byhxhx + byhyhy) - areadivf3 * (dvdx + dvdy) * (dvdx + dvdy);
  H[0][1] = areadivf * (-byhxhx) + areadivf3 * (dvdx + dvdy) * (dvdx);
  H[0][2] = areadivf * (-byhyhy) + areadivf3 * (dvdx + dvdy) * (dvdy);
  H[0][3] = 0.0;
  H[1][1] = areadivf * byhxhx - areadivf3 * dvdx * dvdx;
  H[1][2] = areadivf3 * (-dvdx) * dvdy;
  H[2][2] = areadivf * byhyhy - areadivf3 * dvdy * dvdy;

  /* upper triangle contribution */
  dvdx = (x[3] - x[2]) / hx;
  dvdy = (x[3] - x[1]) / hy;
  fup = sqrt( 1 + dvdx * dvdx + dvdy * dvdy );
  dvdx = dvdx / hx;
  dvdy = dvdy / hy;
  areadivf = area / fup;
  areadivf3 = areadivf / (fup * fup);
  H[1][1] += areadivf * byhyhy - areadivf3 * dvdy * dvdy;
  H[1][2] += areadivf3 * (-dvdy) * dvdx;
  H[2][2] += areadivf * byhxhx - areadivf3 * (dvdx * dvdx);
  H[1][3] = areadivf * (-byhyhy) + areadivf3 * (dvdx + dvdy) * dvdy;
  H[2][3] = areadivf * (-byhxhx) + areadivf3 * (dvdx + dvdy) * dvdx;
  H[3][3] = areadivf * (byhxhx + byhyhy) - areadivf3 * (dvdx + dvdy) * (dvdx + dvdy);

  H[1][0] = H[0][1];
  H[2][0] = H[0][2];
  H[3][0] = H[0][3];
  H[2][1] = H[1][2];
  H[3][1] = H[1][3];
  H[3][2] = H[2][3];

  return 0;

} /* MSurfLocalHessian */


/*------- USER-DEFINED: routine to evaluate the function 
          and gradient at the whole grid             -------*/
#undef __FUNCT__
#define __FUNCT__ "WholeMSurfFunctionGradient"

/*
  WholeMSurfFunctionGradient - Evaluates function and gradient over the 
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
static int WholeMSurfFunctionGradient(TAO_APPLICATION daapplication, DA da, Vec X, double *f, Vec G, void *ptr) {

  AppCtx *user = (AppCtx*)ptr;
  Vec localX, localG;
  PetscInt i, j;
  int info;
  PetscInt xs, xm, gxs, gxm, xe, ys, ym, gys, gym, ye;
  double **x, **g;
  double floc = 0.0;
  PetscScalar zero = 0.0;

  double hx, hy, area;
  double dvdx, dvdy, flow, fup;
  double areadivf;

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
      flow = sqrt( 1 + dvdx * dvdx + dvdy * dvdy );
      areadivf = area / flow;
      g[j][i] += (dvdx / hx + dvdy / hy) * areadivf;
      g[j][i+1] += (-dvdx / hx) * areadivf;
      g[j+1][i] += (-dvdy / hy) * areadivf;

      /* upper triangle contribution */
      dvdx = (x[j+1][i+1] - x[j+1][i]) / hx;
      dvdy = (x[j+1][i+1] - x[j][i+1]) / hy;
      fup = sqrt( 1 + dvdx * dvdx + dvdy * dvdy );
      areadivf = area / fup;
      g[j][i+1] += (-dvdy / hy) * areadivf;
      g[j+1][i] += (-dvdx / hx) * areadivf;
      g[j+1][i+1] += (dvdx / hx + dvdy / hy) * areadivf;

      floc += area * (flow + fup);

    }
  }

  info = MPI_Allreduce(&floc, f, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); CHKERRQ(info);

  info = DAVecRestoreArray(da, localX, (void**)&x); CHKERRQ(info);
  info = DAVecRestoreArray(da, localG, (void**)&g); CHKERRQ(info);

  info = DALocalToGlobalBegin(da, localG, G); CHKERRQ(info);
  info = DALocalToGlobalEnd(da, localG, G); CHKERRQ(info);

  info = DARestoreLocalVector(da, &localX); CHKERRQ(info);
  info = DARestoreLocalVector(da, &localG); CHKERRQ(info);

  info = PetscLogFlops((xe-xs) * (ye-ys) * 42); CHKERRQ(info);

  return 0;
} /* WholeMSurfFunctionGradient  */


/*------- USER-DEFINED: routine to evaluate the Hessian 
          at the whole grid             -------*/

#undef __FUNCT__
#define __FUNCT__ "WholeMSurfHessian"
/*
  WholeMSurfHessian - Evaluates Hessian over the whole grid

  Input:
    daapplication - TAO application object
    da  - distributed array
    X   - the current point, at which the function and gradient are evaluated
    ptr - user-defined application context

  Output:
    H - Hessian at X
*/
static int WholeMSurfHessian(TAO_APPLICATION daapplication, DA da, Vec X, Mat H, void *ptr) {

  AppCtx *user = (AppCtx*)ptr;
  Vec localX;
  int info;
  PetscInt  i, j, ind[4];
  PetscInt xs, xm, gxs, gxm, xe, ys, ym, gys, gym, ye;
  double smallH[4][4];
  double **x;

  double hx, hy, area, byhxhx, byhyhy;
  double dvdx, dvdy, flow, fup;
  double areadivf, areadivf3;
  PetscTruth assembled;

  hx = user->hx;
  hy = user->hy;
  area = user->area;
  
  byhxhx = 1.0 / (hx * hx);
  byhyhy = 1.0 / (hy * hy);

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

      /* 0 is 0,0; 1 is 1,0; 2 is 0,1; 3 is 1,1 */
      dvdx = (x[j][i] - x[j][i+1]) / hx;  /* lower triangle contribution */
      dvdy = (x[j][i] - x[j+1][i]) / hy;
      flow = sqrt( 1 + dvdx * dvdx + dvdy * dvdy );
      dvdx = dvdx / hx;
      dvdy = dvdy / hy;
      areadivf = area / flow;
      areadivf3 = areadivf / (flow * flow);
      smallH[0][0] = areadivf * (byhxhx + byhyhy) - areadivf3 * (dvdx + dvdy) * (dvdx + dvdy);
      smallH[0][1] = areadivf * (-byhxhx) + areadivf3 * (dvdx + dvdy) * (dvdx);
      smallH[0][2] = areadivf * (-byhyhy) + areadivf3 * (dvdx + dvdy) * (dvdy);
      smallH[0][3] = 0.0;
      smallH[1][1] = areadivf * byhxhx - areadivf3 * dvdx * dvdx;
      smallH[1][2] = areadivf3 * (-dvdx) * dvdy;
      smallH[2][2] = areadivf * byhyhy - areadivf3 * dvdy * dvdy;

      /* upper triangle contribution */
      dvdx = (x[j+1][i+1] - x[j+1][i]) / hx;
      dvdy = (x[j+1][i+1] - x[j][i+1]) / hy;
      fup = sqrt( 1 + dvdx * dvdx + dvdy * dvdy );
      dvdx = dvdx / hx;
      dvdy = dvdy / hy;
      areadivf = area / fup;
      areadivf3 = areadivf / (fup * fup);
      smallH[1][1] += areadivf * byhyhy - areadivf3 * dvdy * dvdy;
      smallH[1][2] += areadivf3 * (-dvdy) * dvdx;
      smallH[2][2] += areadivf * byhxhx - areadivf3 * (dvdx * dvdx);
      smallH[1][3] = areadivf * (-byhyhy) + areadivf3 * (dvdx + dvdy) * dvdy;
      smallH[2][3] = areadivf * (-byhxhx) + areadivf3 * (dvdx + dvdy) * dvdx;
      smallH[3][3] = areadivf * (byhxhx + byhyhy) - areadivf3 * (dvdx + dvdy) * (dvdx + dvdy);

      smallH[1][0] = smallH[0][1];
      smallH[2][0] = smallH[0][2];
      smallH[3][0] = smallH[0][3];
      smallH[2][1] = smallH[1][2];
      smallH[3][1] = smallH[1][3];
      smallH[3][2] = smallH[2][3];

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

  info = PetscLogFlops((xe-xs) * (ye-ys) * 83 + 4); CHKERRQ(info);
  return 0;

} /* WholeMSurfHessian */
