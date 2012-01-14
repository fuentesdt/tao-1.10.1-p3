/*$Id: ginzland3.c,v 1.1 2002/05/13 19:32:57 benson Exp $*/

/* Program usage: mpirun -np <proc> ginzland3 [-help] [all TAO options] */

/*
  Include "tao.h" so we can use TAO solvers.
  petscda.h for distributed array
  ad_deriv.h for AD gradient
*/

#include "petscda.h"
#include "tao.h"
#include "taodaapplication.h"


static char  help[] = "";

/*
   User-defined application context - contains data needed by the
   application-provided call-back routines.
*/


/*  
    This structure is used only when an ADIC generated gradient is used.
    An InactiveDouble type is a double 
*/
typedef struct {
  PetscScalar x,y,vpotx,vpoty;
} Field;

typedef struct {
  int                 mx,my;
  InactiveDouble      hx, hy;        /* increment size in both directions */
  InactiveDouble      twopivornum; 
  InactiveDouble      tkappa;
} ADFGCtx;

#ifdef ad_GRAD_MAX
#undef ad_GRAD_MAX
#endif
#define ad_GRAD_MAX 16

int FormFunctionGradient(TAO_APPLICATION,DA,Vec,double *,Vec,void*);
/*
int DAAppSetADElementFunctionGradient2(TAO_APPLICATION, 
				       int (*)(int[2],DERIV_TYPE[16],DERIV_TYPE*,void*), 
				       int, void *);
*/
typedef struct {
  int         fix;
  int         mx,my;
  int         vornum;
  double      tkappa;
  double      hx, hy;        /* increment size in both directions */
  double      twopivornum;
  ADFGCtx     fgctx;         /* Used only when an ADIC generated gradient is used */
} AppCtx;

/* User-defined routines */
static int AppCtxInitialize(void *);
static int AppInitialSolution(TAO_APPLICATION myapp, DA da, Vec X, AppCtx *ctx);
static int MyGridMonitorBefore(TAO_APPLICATION, DA, int, void *);

#undef __FUNCT__
#define __FUNCT__ "main"
int main( int argc, char **argv ) {

  int                  info;         /* used to check for functions returning nonzeros */
  int                  iter,nlevels;
  int                  Nx,Ny;
  double               ff,gnorm;
  DA                   DAarray[20];
  PetscTruth           flg;                             /* flags */
  TaoMethod            method = "tao_nls";                      /* minimization method */
  TaoTerminateReason   reason;
  TAO_SOLVER           tao;                               /* TAO_SOLVER solver context */
  TAO_APPLICATION      GinzburghLandau;                       /* The PETSc application */
  AppCtx               user;                              /* user-defined work context */
  Vec x;
  double zero=0.0;
  KSP    ksp;
  PC     pc;

  /* Initialize TAO */
  PetscInitialize(&argc, &argv, (char *)0, help);
  TaoInitialize(&argc, &argv, (char *)0, help);

  info = AppCtxInitialize((void*)&user); CHKERRQ(info);

  nlevels=1;
  info = PetscOptionsGetInt(PETSC_NULL,"-nlevels",&nlevels,&flg); CHKERRQ(info);

  PetscPrintf(MPI_COMM_WORLD,"\n---- Ginzburg - Landau Superconductivity -----\n\n");

  /* Let PETSc determine the vector distribution */
  Nx = PETSC_DECIDE; Ny = PETSC_DECIDE;

  /* Create distributed array (DA) to manage parallel grid and vectors  */
  info = DACreate2d(PETSC_COMM_WORLD,DA_XYPERIODIC, DA_STENCIL_BOX, user.mx,
                    user.my,Nx,Ny,4,1,PETSC_NULL,PETSC_NULL,&DAarray[0]); CHKERRQ(info);
  for (iter=1;iter<nlevels;iter++){
    info = DARefine(DAarray[iter-1],PETSC_COMM_WORLD,&DAarray[iter]); CHKERRQ(info);
  }

  /* Create TAO solver and set desired solution method */
  info = TaoCreate(MPI_COMM_WORLD,method,&tao); CHKERRQ(info);
  info = TaoApplicationCreate(PETSC_COMM_WORLD, &GinzburghLandau); CHKERRQ(info);
  info = TaoAppSetDAApp(GinzburghLandau,DAarray,nlevels); CHKERRQ(info);
  info = DAAppSetBeforeMonitor(GinzburghLandau,MyGridMonitorBefore,(void*)&user); CHKERRQ(info);
  info = DAAppPrintInterpolationError(GinzburghLandau); CHKERRQ(info);
  info = DAAppPrintStageTimes(GinzburghLandau); CHKERRQ(info);

  info = TaoAppSetRelativeTolerance(GinzburghLandau,1.0e-6); CHKERRQ(info);
  info = TaoSetTolerances(tao,0,0,0,0); CHKERRQ(info);
  info = TaoSetGradientTolerances(tao,0,0,0); CHKERRQ(info);
  info = DAAppSetObjectiveAndGradientRoutine(GinzburghLandau,FormFunctionGradient,(void*)&user); CHKERRQ(info);

  //  info = DAAppSetElementObjectiveAndGradientRoutine(GinzburghLandau,GinzLandLocalFunctionGradient,0,(void*)&user); CHKERRQ(info);
  //  info = DAAppSetADElementFunctionGradient2(GinzburghLandau,ad_GinzLandLocalFunction,150,(void *)&user.fgctx); CHKERRQ(info);

  info = DAAppSetHessianMat(GinzburghLandau); CHKERRQ(info);
  info = TaoAppSetHessianRoutine(GinzburghLandau,TaoAppDefaultComputeHessianColor,(void*)GinzburghLandau);CHKERRQ(info);

  /* Check for any tao command line options */
  info = DAAppGetSolution(GinzburghLandau,0,&x); CHKERRQ(info);
  info = VecSet(x, zero); CHKERRQ(info);
  info = AppInitialSolution(GinzburghLandau,DAarray[0], x, &user); CHKERRQ(info);
  //  info=VecView(x,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(info);
  //  info = DAAppSetInitialSolution(GinzburghLandau, x); CHKERRQ(info);

  info = TaoAppGetKSP(GinzburghLandau,&ksp);CHKERRQ(info);
  info = KSPSetType(ksp,KSPCG); CHKERRQ(info);
  info = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,1000);CHKERRQ(info);
  info = KSPGetPC(ksp,&pc);CHKERRQ(info);
  info = PCSetType(pc,PCJACOBI);CHKERRQ(info);

  info = TaoSetOptions(GinzburghLandau,tao); CHKERRQ(info);

  /* SOLVE THE APPLICATION */
  info = TaoDAAppSolve(GinzburghLandau,tao);  CHKERRQ(info);

  /* Get information on termination */
  info = TaoGetSolutionStatus(tao,&iter,&ff,&gnorm,0,0,&reason); CHKERRQ(info);
  if (reason <= 0 ){
    PetscPrintf(MPI_COMM_WORLD,"Try a different TAO method, adjust some parameters, or check the function evaluation routines\n");
    PetscPrintf(MPI_COMM_WORLD," Iterations: %d,  Function Value: %4.2e, Residual: %4.2e \n",iter,ff,gnorm);
  }

  info = PetscOptionsHasName(PETSC_NULL,"-view_sol",&flg); CHKERRQ(info);
  if (flg){
    info = DAAppGetSolution(GinzburghLandau,nlevels-1,&x); CHKERRQ(info);
    info = VecView(x,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(info);
  }

  /*  To View TAO solver information */
  info = TaoView(tao); CHKERRQ(info);

  /* Free TAO data structures */
  info = TaoDestroy(tao); CHKERRQ(info);
  info = TaoAppDestroy(GinzburghLandau); CHKERRQ(info);

  /* Free PETSc data structures */
  for (iter=0;iter<nlevels;iter++){
    info = DADestroy(DAarray[iter]); CHKERRQ(info);
  }

  /* Finalize TAO */
  TaoFinalize();
  PetscFinalize();

  return 0;
} /* main */

static int AppInitialSolution(TAO_APPLICATION myapp, DA da, Vec X, AppCtx *ctx) {

  int xs, xm, mx, ys, ym, my;
  double tkappa, hx, hy, sqn, ppi,bave, sqrtv, vornum;
  double xpt,ypt,t1;
  Field **x;
  int i,j, info;

  vornum=ctx->vornum;
  tkappa=ctx->tkappa;

  info = DAGetInfo(da,PETSC_NULL,&mx,&my,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(info);
  info = DAGetCorners(da, &xs, &ys, PETSC_NULL, &xm, &ym, PETSC_NULL); CHKERRQ(info);

  info = MyGridMonitorBefore(myapp,da,0,(void*)ctx); CHKERRQ(info);

  tkappa=ctx->tkappa; vornum=ctx->vornum;
  hx=ctx->hx; hy=ctx->hy;
  sqn=mx*my; ppi=4.0*atan(1.0);
  bave=2.0*(ppi)*(vornum)*(tkappa)/(sqn*hx*hy);
  sqrtv=sqrt(vornum*1.0)*(ppi);

  info = DAVecGetArray(da, X, (void**)&x); CHKERRQ(info);
  for (j = ys; j < ys+ym; j++) {
    ypt=j*hy;
    for (i = xs; i < xs+xm; i++) {
      xpt=i*hx;
      t1= sin(sqrtv*xpt/6.0) * sin(sqrtv*ypt/(6.0*sqrt(3.0))) ;
      x[j][i].x = 1.0 - t1*t1;
      x[j][i].y = 0.0;
    }
  }

  for (j = ys; j < ys+ym; j++) {
    for (i = xs; i < xs+xm; i++) {
      xpt=i*hx;
      x[j][i].vpotx = 0.0;
      x[j][i].vpoty = bave*xpt/TaoMax(tkappa, 1e-4);
    }
  }
  info = DAVecRestoreArray(da, X, (void**)&x); CHKERRQ(info);

  return 0;
}

/*----- The following two routines
  MyGridMonitorBefore    MyGridMonitorAfter
  help diplay info of iterations at every grid level 
*/

#undef __FUNCT__
#define __FUNCT__ "MyGridMonitorBefore"
static int MyGridMonitorBefore(TAO_APPLICATION myapp, DA da, int level, void *ctx) {

  AppCtx *user = (AppCtx*)ctx;
  int info,mx,my;
  double vornum;
  KSP ksp;

  info = DAGetInfo(da,PETSC_NULL,&mx,&my,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,
		   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(info);

  vornum=user->vornum;
  user->hx=sqrt(vornum/2.0)*3.0/(mx);
  user->hy=sqrt(vornum/2.0)*3.0*sqrt(3.0)/(my);
  user->twopivornum = 2.0*(4.0*atan(1.0))*vornum;
  user->mx = mx;
  user->my = my;

  user->fgctx.hx          = user->hx;
  user->fgctx.hy          = user->hy;
  user->fgctx.mx          = user->mx;
  user->fgctx.my          = user->my;
  user->fgctx.tkappa      = user->tkappa;
  user->fgctx.twopivornum = user->twopivornum;

  info = TaoAppGetKSP(myapp,&ksp);CHKERRQ(info);
  if (mx<50){
    info = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,500);CHKERRQ(info);
  } else {
    info = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,1000);CHKERRQ(info);
  }

  if (level > 0) {
    user->fix = 1;
  }
  else {
    user->fix = 0;
  }
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
  user->mx = user->my = 12;
  user->vornum=8;
  user->tkappa=5.0;

  /* Check for command line arguments that override defaults */
  info = PetscOptionsGetInt(TAO_NULL, "-mx", &user->mx, &flg); CHKERRQ(info);
  info = PetscOptionsGetInt(TAO_NULL, "-my", &user->my, &flg); CHKERRQ(info);
  info = PetscOptionsGetInt(TAO_NULL, "-par", &user->vornum, &flg); CHKERRQ(info);
  info = PetscOptionsGetReal(TAO_NULL, "-kappa", &user->tkappa, &flg); CHKERRQ(info);

  user->hx=sqrt(user->vornum/2.0)*3.0/(user->mx);
  user->hy=sqrt(user->vornum/2.0)*3.0*sqrt(3.0)/(user->my);

  user->fgctx.hx          = user->hx;
  user->fgctx.hy          = user->hy;
  user->fgctx.mx          = user->mx;
  user->fgctx.my          = user->my;
  user->fgctx.tkappa      = user->tkappa;
  user->fgctx.twopivornum = 2.0*(4.0*atan(1.0))*(user->vornum);

  PetscLogFlops(8);

  return 0;

} /* AppCtxInitialize */

#undef __FUNCT__
#define __FUNCT__ "FormFunctionGradient"
int FormFunctionGradient(TAO_APPLICATION app, DA da, Vec X, double *f, Vec G, void *ctx)
{
  AppCtx * user = (AppCtx *) ctx;
  int    info,i,j;
  double floc=0,fcond,fkin,ffield;
  double arg,cosarg,sinarg;
  int    mx=user->mx, my=user->my;
  double hx=user->hx, hy=user->hy;
  double hxhx = hx*hx;
  double hyhy = hy*hy;
  int    xs,xm,gxs,gxm,ys,ym,gys,gym;
  double tkappa=user->tkappa;
  double twopivornum=user->twopivornum;
  double d,d1,d2,delsq,sqn=mx*my,twooversqn=2.0/sqn;
  Field **g;
  Field **x;
  double x4, x5, x7,xy;
  double c1,c2,e1,e2,cosc1,sinc1,cosc2,sinc2;
  Vec    localX,localG;

  PetscFunctionBegin;
  /* Get local mesh boundaries */
  info = DAGetCorners(da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL); CHKERRQ(info);
  info = DAGetGhostCorners(da,&gxs,&gys,PETSC_NULL,&gxm,&gym,PETSC_NULL); CHKERRQ(info);

  /* Scatter ghost points to local vector */
  info = DAGetLocalVector(da, &localX); CHKERRQ(info);
  info = DAGetLocalVector(da, &localG); CHKERRQ(info);
  info = VecSet(G, 0.0); CHKERRQ(info);
  info = VecSet(localG, 0.0); CHKERRQ(info);

  info = DAGlobalToLocalBegin(da, X, INSERT_VALUES, localX); CHKERRQ(info);
  info = DAGlobalToLocalEnd(da, X, INSERT_VALUES, localX); CHKERRQ(info);

  /* Get pointers to vector data */
  info = DAVecGetArray(da,localX,(void**)&x);
  info = DAVecGetArray(da,localG,(void**)&g);

  if (user->fix) {
    // Fix the defect
    double xv, yv;

    // Update boundary conditions for interpolation
    for (j = ys; j <= ys+ym; ++j) {
      for (i = xs; i < xs+xm; ++i) {
        if (i == mx-1) {
	  arg = (j*twopivornum/my);
	  cosarg = cos(arg); 
          sinarg = sin(arg);
  
          xv = x[j][i+1].x;
          yv = x[j][i+1].y;
  
	  x[j][i+1].x = xv*cosarg - yv*sinarg;
	  x[j][i+1].y = xv*sinarg + yv*cosarg;
	  x[j][i+1].vpotx = x[j][i+1].vpotx;
	  x[j][i+1].vpoty = x[j][i+1].vpoty + (twopivornum) / (my*hy);
        }
      }
    }
 
    for (j = ys; j < ys+ym; ++j) {
      if (0 == j % 2) {
        for (i = xs; i < xs+xm; ++i) {
          if (0 == i % 2) {
            // Coarse point -- no interpolation
          }
          else {
            // Fine point -- x interpolation
            x[j][i].x = (x[j][i-1].x + x[j][i+1].x) / 2.0;
            x[j][i].y = (x[j][i-1].y + x[j][i+1].y) / 2.0;
            x[j][i].vpotx = (x[j][i-1].vpotx + x[j][i+1].vpotx) / 2.0;
            x[j][i].vpoty = (x[j][i-1].vpoty + x[j][i+1].vpoty) / 2.0;
          }
        }
      }
      else {
        for (i = xs; i < xs+xm; ++i) {
          if (0 == i % 2) {
            // Fine point -- y interpolation
            x[j][i].x = (x[j-1][i].x + x[j+1][i].x) / 2.0;
            x[j][i].y = (x[j-1][i].y + x[j+1][i].y) / 2.0;
            x[j][i].vpotx = (x[j-1][i].vpotx + x[j+1][i].vpotx) / 2.0;
            x[j][i].vpoty = (x[j-1][i].vpoty + x[j+1][i].vpoty) / 2.0;
          }
          else {
            // Fine point -- x,y interpolation
            x[j][i].x = (x[j-1][i-1].x + x[j-1][i+1].x + x[j+1][i-1].x + x[j+1][i+1].x) / 4.0;
            x[j][i].y = (x[j-1][i-1].y + x[j-1][i+1].y + x[j+1][i-1].y + x[j+1][i+1].y) / 4.0;
            x[j][i].vpotx = (x[j-1][i-1].vpotx + x[j-1][i+1].vpotx + x[j+1][i-1].vpotx + x[j+1][i+1].vpotx) / 4.0;
            x[j][i].vpoty = (x[j-1][i-1].vpoty + x[j-1][i+1].vpoty + x[j+1][i-1].vpoty + x[j+1][i+1].vpoty) / 4.0;
          }
        }
      }
    }

#if 1
    // Coarse function evaluation

    double c_fcond = 0;
    double c_fkin1 = 0;
    double c_fkin2 = 0;
    double c_ffield = 0;

    for (j=ys; j<ys+ym; j++) {
      if (0 == j % 2) {
        for (i=xs; i<xs+xm;i++) {
          if (0 == i % 2) {
            delsq = x[j][i].x*x[j][i].x + x[j][i].y*x[j][i].y;
            c_fcond += (-delsq + (delsq*delsq)/2.0 );

            c1 = 2.0*hx*x[j][i].vpotx;
            cosc1 = cos(c1);
            sinc1 = sin(c1);

            d1 = x[j][i+2].x - x[j][i].x*cosc1 + x[j][i].y*sinc1;
            d2 = x[j][i+2].y - x[j][i].y*cosc1 - x[j][i].x*sinc1;

            c2 = 2.0*hy*x[j][i].vpoty;
            cosc2 = cos(c2);
            sinc2 = sin(c2);

            e1 = x[j+2][i].x - x[j][i].x*cosc2 + x[j][i].y*sinc2;
            e2 = x[j+2][i].y - x[j][i].y*cosc2 - x[j][i].x*sinc2;

            c_fkin1 += (d1*d1 + d2*d2)/(4.0*hxhx);
            c_fkin2 += (e1*e1 + e2*e2)/(4.0*hyhy);

            // d1 = (x[j][i+2].x - x[j][i].x) / (2.0 * hx) + x[j][i].vpotx * x[j][i].y;
            // d2 = (x[j][i+2].y - x[j][i].y) / (2.0 * hx) - x[j][i].vpotx * x[j][i].x;
	
            // e1 = (x[j+2][i].x - x[j][i].x) / (2.0 * hy) + x[j][i].vpoty * x[j][i].y;
            // e2 = (x[j+2][i].y - x[j][i].y) / (2.0 * hy) - x[j][i].vpoty * x[j][i].x;

            // c_fkin1 += d1*d1 + d2*d2;
            // c_fkin2 += e1*e1 + e2*e2;

            d1 = (x[j][i+2].vpoty - x[j][i].vpoty)/(2.0*hx) - (x[j+2][i].vpotx - x[j][i].vpotx)/(2.0*hy);
            c_ffield += tkappa * tkappa * d1 * d1;
          }
        }
      }
    }

    // Fine function evaluation

    double f_fcond = 0;
    double f_fkin1 = 0;
    double f_fkin2 = 0;
    double f_ffield = 0;

    for (j=ys; j<ys+ym; j++) {
      for (i=xs; i<xs+xm;i++) {
        delsq = x[j][i].x*x[j][i].x + x[j][i].y*x[j][i].y;
        f_fcond += (-delsq + (delsq*delsq)/2.0 );

        c1 = hx*x[j][i].vpotx;
        cosc1 = cos(c1);
        sinc1 = sin(c1);

        d1 = x[j][i+1].x - x[j][i].x*cosc1 + x[j][i].y*sinc1;
        d2 = x[j][i+1].y - x[j][i].y*cosc1 - x[j][i].x*sinc1;

        c2 = hy*x[j][i].vpoty;
        cosc2 = cos(c2);
        sinc2 = sin(c2);

        e1 = x[j+1][i].x - x[j][i].x*cosc2 + x[j][i].y*sinc2;
        e2 = x[j+1][i].y - x[j][i].y*cosc2 - x[j][i].x*sinc2;

        f_fkin1 += (d1*d1 + d2*d2)/(hxhx);
        f_fkin2 += (e1*e1 + e2*e2)/(hyhy);

        // d1 = (x[j][i+1].x - x[j][i].x) / hx + x[j][i].vpotx * x[j][i].y;
        // d2 = (x[j][i+1].y - x[j][i].y) / hx - x[j][i].vpotx * x[j][i].x;

        // e1 = (x[j+1][i].x - x[j][i].x) / hy + x[j][i].vpoty * x[j][i].y;
        // e2 = (x[j+1][i].y - x[j][i].y) / hy - x[j][i].vpoty * x[j][i].x;

        // f_fkin1 += d1*d1 + d2*d2;
        // f_fkin2 += e1*e1 + e2*e2;

        xy = (x[j][i+1].vpoty - x[j][i].vpoty)/(hx) - (x[j+1][i].vpotx - x[j][i].vpotx)/(hy);
        d = tkappa*tkappa*xy;
        f_ffield += d*xy;
      }
    }

    c_fcond = 4.0 * c_fcond / sqn;
    c_fkin1 = 4.0 * c_fkin1 / sqn;
    c_fkin2 = 4.0 * c_fkin2 / sqn;
    c_ffield = 4.0 * c_ffield / sqn;

    f_fcond = f_fcond / sqn;
    f_fkin1 = f_fkin1 / sqn;
    f_fkin2 = f_fkin2 / sqn;
    f_ffield = f_ffield / sqn;

    printf("%5.4e %5.4e %5.4e %5.4e\n",
           fabs(c_fcond - f_fcond),
           fabs(c_fkin1 - f_fkin1),
           fabs(c_fkin2 - f_fkin2),
           fabs(c_ffield - f_ffield));
           
    printf("%5.4e %5.4e\n", 
           c_fcond + c_fkin1 + c_fkin2 + c_ffield,
           f_fcond + f_fkin1 + f_fkin2 + f_ffield);
#endif

    info = DAVecRestoreArray(da, localX, (void**)&x); CHKERRQ(info);
    info = DALocalToGlobal(da, localX, INSERT_VALUES, X); CHKERRQ(info);
    info = DARestoreLocalVector(da, &localX); CHKERRQ(info);

    info = DAGetLocalVector(da, &localX); CHKERRQ(info);
    info = DAGlobalToLocalBegin(da, X, INSERT_VALUES, localX); CHKERRQ(info);
    info = DAGlobalToLocalEnd(da, X, INSERT_VALUES, localX); CHKERRQ(info);
    info = DAVecGetArray(da, localX, (void**)&x);
    user->fix = 0;
  }

  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm;i++) {
      if (i == mx-1) {
	arg = (j*twopivornum/my);
	cosarg = cos(arg); sinarg = sin(arg);
	x4 = x[j][i+1].x*cosarg - x[j][i+1].y*sinarg;
	x5 = x[j][i+1].x*sinarg + x[j][i+1].y*cosarg;
	x7 = x[j][i+1].vpoty + (twopivornum) / (my*hy);
      }
      else {
	x4=x[j][i+1].x;
	x5=x[j][i+1].y;
	x7=x[j][i+1].vpoty;
      }

      /*  Compute the Condensation Energy Density*/
      delsq = x[j][i].x * x[j][i].x + x[j][i].y * x[j][i].y;
      fcond = delsq * delsq / 2.0 - delsq;

      g[j][i].x += twooversqn * x[j][i].x * (delsq - 1.0);
      g[j][i].y += twooversqn * x[j][i].y * (delsq - 1.0);

      /*  Compute the Kinetic Energy Density. */
      c1 = hx*x[j][i].vpotx;
      cosc1 = cos(c1);
      sinc1 = sin(c1);
      d1 = x4 - x[j][i].x*cosc1 + x[j][i].y*sinc1;
      d2 = x5 - x[j][i].y*cosc1 - x[j][i].x*sinc1;

      c2 = hy*x[j][i].vpoty;
      cosc2 = cos(c2);
      sinc2 = sin(c2);
      e1 = x[j+1][i].x - x[j][i].x*cosc2 + x[j][i].y*sinc2;
      e2 = x[j+1][i].y - x[j][i].y*cosc2 - x[j][i].x*sinc2;

      fkin = (d1*d1 + d2*d2)/(hxhx) + (e1*e1 + e2*e2)/(hyhy);

      // g[j][i].x -= twooversqn * (d1*cosc1+d2*sinc1)/(hxhx);
      // g[j][i].x -= twooversqn * (e1*cosc2+e2*sinc2)/(hyhy);
      // g[j][i].y += twooversqn * (d1*sinc1-d2*cosc1)/(hxhx);
      // g[j][i].y += twooversqn * (e1*sinc2-e2*cosc2)/(hyhy);

      // g[j][i].vpotx += twooversqn / hx * 
      //   (d1*(x[j][i].x * sinc1 + x[j][i].y*cosc1) +
      //    d2*(x[j][i].y * sinc1 - x[j][i].x * cosc1));

      // g[j][i].vpoty += twooversqn / hy * 
      //   (e1*(x[j][i].x * sinc2 + x[j][i].y * cosc2) +
      //    e2*(x[j][i].y * sinc2 - x[j][i].x * cosc2));

      // if (i == mx-1) {
      //   g[j][i+1].x += twooversqn / hxhx * (d1*cosarg + d2*sinarg);
      //   g[j][i+1].y += twooversqn / hxhx * (d2*cosarg - d1*sinarg);
      // }
      // else {
      //   g[j][i+1].x += twooversqn / hxhx * d1;
      //   g[j][i+1].y += twooversqn / hxhx * d2;
      // }

      // g[j+1][i].x += twooversqn / hyhy * e1;
      // g[j+1][i].y += twooversqn / hyhy * e2;

      d1 = (x4 - x[j][i].x) / hx + x[j][i].vpotx * x[j][i].y;
      d2 = (x5 - x[j][i].y) / hx - x[j][i].vpotx * x[j][i].x;

      e1 = (x[j+1][i].x - x[j][i].x) / hy + x[j][i].vpoty * x[j][i].y;
      e2 = (x[j+1][i].y - x[j][i].y) / hy - x[j][i].vpoty * x[j][i].x;

      fkin = d1*d1 + d2*d2 + e1*e1 + e2*e2;

      g[j][i].x -= twooversqn * (d1 / hx + d2 * x[j][i].vpotx + 
                                 e1 / hy + e2 * x[j][i].vpoty);
 
      g[j][i].y -= twooversqn * (d2 / hx - d1 * x[j][i].vpotx +
                                 e2 / hy - e1 * x[j][i].vpoty);
      g[j][i].vpotx += twooversqn * (d1 * x[j][i].y - d2 * x[j][i].x);
      g[j][i].vpoty += twooversqn * (e1 * x[j][i].y - e2 * x[j][i].x);

      if (i == mx - 1) {
        g[j][i+1].x += twooversqn * (d1 * cosarg + d2 * sinarg) / hx;
        g[j][i+1].y += twooversqn * (d2 * cosarg - d1 * sinarg) / hx;
      }
      else {
        g[j][i+1].x += twooversqn * d1 / hx;
        g[j][i+1].y += twooversqn * d2 / hx;
      }

      g[j+1][i].x += twooversqn * e1 / hy;
      g[j+1][i].y += twooversqn * e2 / hy;

      /*  Compute the Magnetic Field Energy Density. */
      d1 = (x7 - x[j][i].vpoty) / hx - (x[j+1][i].vpotx - x[j][i].vpotx) / hy;
      ffield = tkappa * tkappa * d1 * d1;

      g[j][i].vpotx += twooversqn * tkappa * tkappa * d1 / hy;
      g[j][i].vpoty -= twooversqn * tkappa * tkappa * d1 / hx;

      g[j+1][i].vpotx -= twooversqn * tkappa * tkappa * d1 / hy;
      g[j][i+1].vpoty += twooversqn * tkappa * tkappa * d1 / hx;

      floc += (fcond + fkin + ffield)/sqn;
    }
  }

  info = MPI_Allreduce(&floc, f, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); CHKERRQ(info);

  info = DAVecRestoreArray(da, localX, (void**)&x); CHKERRQ(info);
  info = DAVecRestoreArray(da, localG, (void**)&g); CHKERRQ(info);
      
  info = DALocalToGlobalBegin(da, localG, G); CHKERRQ(info);
  info = DALocalToGlobalEnd(da, localG, G); CHKERRQ(info);

  info = DARestoreLocalVector(da, &localX); CHKERRQ(info);
  info = DARestoreLocalVector(da, &localG); CHKERRQ(info);
  PetscFunctionReturn(0);
}
