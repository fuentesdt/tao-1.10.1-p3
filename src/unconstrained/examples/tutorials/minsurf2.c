/*$Id$*/

/* Program usage: mpirun -np <proc> minsurf2 [-help] [all TAO options] */

/*
  Include "tao.h" so we can use TAO solvers.
  petscda.h for distributed array
*/
#include "petscda.h"
#include "tao.h"

static  char help[] = 
"This example demonstrates use of the TAO package to \n\
solve an unconstrained minimization problem.  This example is based on a \n\
problem from the MINPACK-2 test suite.  Given a rectangular 2-D domain and \n\
boundary values along the edges of the domain, the objective is to find the\n\
surface with the minimal area that satisfies the boundary conditions.\n\
The command line options are:\n\
  -mx <xg>, where <xg> = number of grid points in the 1st coordinate direction\n\
  -my <yg>, where <yg> = number of grid points in the 2nd coordinate direction\n\
  -start <st>, where <st> =0 for zero vector, <st> >0 for random start, and <st> <0 \n\
               for an average of the boundary conditions\n\n";

/*T
   Concepts: TAO - Solving an unconstrained minimization problem
   Routines: TaoInitialize(); TaoFinalize();
   Routines: TaoCreate(); TaoDestroy();
   Routines: TaoApplicationCreate(); TaoAppDestroy();
   Routines: TaoAppSetInitialSolutionVec();
   Routines: TaoAppSetObjectiveAndGradientRoutine();
   Routines: TaoAppSetHessianMat(); TaoAppSetHessianRoutine();
   Routines: TaoSetOptions();
   Routines: TaoAppGetKSP(); TaoSolveApplication();
   Routines: TaoAppSetMonitor(); TaoView();
   Routines: TaoAppGetSolutionVec();
   Processors: 1
T*/

/* 
   User-defined application context - contains data needed by the 
   application-provided call-back routines, FormFunctionGradient() 
   and FormHessian().
*/
typedef struct {
  PetscInt      mx, my;                 /* discretization in x, y directions */
  double      *bottom, *top, *left, *right;             /* boundary values */
  DA          da;                      /* distributed array data structure */
  Mat         H;                       /* Hessian */
  ISColoring  iscoloring;
} AppCtx;


/* -------- User-defined Routines --------- */

static int MSA_BoundaryConditions(AppCtx*);
static int MSA_InitialPoint(AppCtx*,Vec);
int QuadraticH(AppCtx*,Vec,Mat);
int FormFunctionGradient(TAO_APPLICATION,Vec,double *,Vec,void*);
int FormGradient(TAO_APPLICATION,Vec,Vec,void*);
int FormHessian(TAO_APPLICATION,Vec,Mat*,Mat*,MatStructure *,void*);
int My_Monitor(TAO_APPLICATION, void *);

#undef __FUNCT__
#define __FUNCT__ "main"
int main( int argc, char **argv )
{
  int             info;                /* used to check for functions returning nonzeros */
  PetscInt          Nx, Ny;              /* number of processors in x- and y- directions */
  PetscInt          iter;                /* iteration information */
  double          ff,gnorm;
  Vec             x;                   /* solution, gradient vectors */
  PetscTruth      flg, viewmat;        /* flags */
  PetscTruth      fddefault, fdcoloring;   /* flags */
  KSP             ksp;                 /* Krylov subspace method */
  TaoMethod       method = "tao_nls";  /* minimization method */
  TaoTerminateReason reason;           
  TAO_SOLVER      tao;                 /* TAO_SOLVER solver context */
  TAO_APPLICATION minsurfapp;          /* The PETSc application */
  AppCtx          user;                /* user-defined work context */

  /* Initialize TAO */
  PetscInitialize( &argc, &argv,(char *)0,help );
  TaoInitialize( &argc, &argv,(char *)0,help );

  /* Specify dimension of the problem */
  user.mx = 10; user.my = 10;

  /* Check for any command line arguments that override defaults */
  info = PetscOptionsGetInt(PETSC_NULL,"-mx",&user.mx,&flg); CHKERRQ(info);
  info = PetscOptionsGetInt(PETSC_NULL,"-my",&user.my,&flg); CHKERRQ(info);

  PetscPrintf(MPI_COMM_WORLD,"\n---- Minimum Surface Area Problem -----\n");
  PetscPrintf(MPI_COMM_WORLD,"mx: %d     my: %d   \n\n",user.mx,user.my);


  /* Let PETSc determine the vector distribution */
  Nx = PETSC_DECIDE; Ny = PETSC_DECIDE;

  /* Create distributed array (DA) to manage parallel grid and vectors  */
  info = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,user.mx,
                    user.my,Nx,Ny,1,1,PETSC_NULL,PETSC_NULL,&user.da); CHKERRQ(info);
  

  /* Create TAO solver and set desired solution method. Create an TAO application structure */
  info = TaoCreate(PETSC_COMM_WORLD,method,&tao); CHKERRQ(info);
  info = TaoApplicationCreate(PETSC_COMM_WORLD,&minsurfapp); CHKERRQ(info);

  /*
     Extract global vector from DA for the vector of variables --  PETSC routine
     Compute the initial solution                              --  application specific, see below
     Set this vector for use by TAO                            --  TAO routine
  */
  info = DACreateGlobalVector(user.da,&x); CHKERRQ(info);
  info = MSA_BoundaryConditions(&user); CHKERRQ(info);         
  info = MSA_InitialPoint(&user,x); CHKERRQ(info);
  info = TaoAppSetInitialSolutionVec(minsurfapp,x); CHKERRQ(info);

  /* 
     Initialize the Application context for use in function evaluations  --  application specific, see below.
     Set routines for function and gradient evaluation 
  */
  info = TaoAppSetObjectiveAndGradientRoutine(minsurfapp,FormFunctionGradient,(void *)&user); CHKERRQ(info);

  /* 
     Given the command line arguments, calculate the hessian with either the user-
     provided function FormHessian, or the default finite-difference driven Hessian
     functions 
  */
  info = PetscOptionsHasName(PETSC_NULL,"-tao_fddefault",&fddefault);CHKERRQ(info);
  info = PetscOptionsHasName(PETSC_NULL,"-tao_fdcoloring",&fdcoloring);CHKERRQ(info);


  /* 
     Create a matrix data structure to store the Hessian and set 
     the Hessian evalution routine.
     Set the matrix structure to be used for Hessian evalutions
  */
  info = DAGetMatrix(user.da,MATAIJ,&user.H);CHKERRQ(info);
  info = MatSetOption(user.H,MAT_SYMMETRIC,PETSC_TRUE); CHKERRQ(info);

  info = TaoAppSetHessianMat(minsurfapp,user.H,user.H); CHKERRQ(info);

  if (fdcoloring) {
    info = DAGetColoring(user.da,IS_COLORING_GLOBAL,MATAIJ,&(user.iscoloring));
    CHKERRQ(info);
    info = TaoAppSetColoring(minsurfapp, user.iscoloring); CHKERRQ(info);
    info = TaoAppSetHessianRoutine(minsurfapp,TaoAppDefaultComputeHessianColor,(void *)&user); CHKERRQ(info);
  } else if (fddefault){
    info = TaoAppSetHessianRoutine(minsurfapp,TaoAppDefaultComputeHessian,(void *)&user); CHKERRQ(info);
  } else {
    info = TaoAppSetHessianRoutine(minsurfapp,FormHessian,(void *)&user); CHKERRQ(info);
  }


  /* 
     If my_monitor option is in command line, then use the user-provided
     monitoring function
  */
  info = PetscOptionsHasName(PETSC_NULL,"-my_monitor",&viewmat); CHKERRQ(info);
  if (viewmat){
    info = TaoAppSetMonitor(minsurfapp,My_Monitor,TAO_NULL); CHKERRQ(info);
  }

  /* Check for any tao command line options */
  info = TaoSetOptions(minsurfapp,tao); CHKERRQ(info);

  /* Limit the number of iterations in the KSP linear solver */
  info = TaoAppGetKSP(minsurfapp,&ksp); CHKERRQ(info);
  if (ksp) {                                              /* Modify the PETSc KSP structure */
    info = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,user.mx*user.my);
    CHKERRQ(info);
  }

  /* SOLVE THE APPLICATION */
  info = TaoSolveApplication(minsurfapp,tao);  CHKERRQ(info);

  /* Get information on termination */
  info = TaoGetSolutionStatus(tao,&iter,&ff,&gnorm,0,0,&reason); CHKERRQ(info);
  if (reason <= 0 ){
    PetscPrintf(MPI_COMM_WORLD,"Try a different TAO method, adjust some parameters, or check the function evaluation routines\n");
    PetscPrintf(MPI_COMM_WORLD," Iterations: %d,  Function Value: %4.2e, Residual: %4.2e \n",iter,ff,gnorm);
  }

  /* 
     To View TAO solver information use
      info = TaoView(tao); CHKERRQ(info); 
  */

  /* Free TAO data structures */
  info = TaoDestroy(tao); CHKERRQ(info);
  info = TaoAppDestroy(minsurfapp); CHKERRQ(info);

  /* Free PETSc data structures */
  info = VecDestroy(x); CHKERRQ(info);
  info = MatDestroy(user.H); CHKERRQ(info);
  if (fdcoloring) {
    info = ISColoringDestroy(user.iscoloring);CHKERRQ(info);
  }
  info = PetscFree(user.bottom); CHKERRQ(info);
  info = PetscFree(user.top); CHKERRQ(info);
  info = PetscFree(user.left); CHKERRQ(info);
  info = PetscFree(user.right); CHKERRQ(info);
  info = DADestroy(user.da); CHKERRQ(info);

  /* Finalize TAO */
  TaoFinalize();
  PetscFinalize();
  
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "FormGradient"
int FormGradient(TAO_APPLICATION taoapp, Vec X, Vec G,void *userCtx){
  int info;
  double fcn;
  TaoFunctionBegin;
  info = FormFunctionGradient(taoapp,X,&fcn,G,userCtx);CHKERRQ(info);
  TaoFunctionReturn(0);
}

/* -------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormFunctionGradient"
/*  FormFunctionGradient - Evaluates the function and corresponding gradient.

    Input Parameters:
.   taoapp     - the TAO_APPLICATION context
.   XX      - input vector
.   userCtx - optional user-defined context, as set by TaoSetFunctionGradient()
    
    Output Parameters:
.   fcn     - the newly evaluated function
.   GG       - vector containing the newly evaluated gradient
*/
int FormFunctionGradient(TAO_APPLICATION taoapp, Vec X, double *fcn,Vec G,void *userCtx){

  AppCtx * user = (AppCtx *) userCtx;
  int    info;
  PetscInt i,j;
  PetscInt mx=user->mx, my=user->my;
  PetscInt xs,xm,gxs,gxm,ys,ym,gys,gym;
  double ft=0;
  double hx=1.0/(mx+1),hy=1.0/(my+1), hydhx=hy/hx, hxdhy=hx/hy, area=0.5*hx*hy;
  double rhx=mx+1, rhy=my+1;
  double f1,f2,f3,f4,f5,f6,d1,d2,d3,d4,d5,d6,d7,d8,xc,xl,xr,xt,xb,xlt,xrb;
  double df1dxc,df2dxc,df3dxc,df4dxc,df5dxc,df6dxc;
  PetscScalar **g, **x;
  Vec    localX;

  /* Get local mesh boundaries */
  info = DAGetLocalVector(user->da,&localX);CHKERRQ(info);

  info = DAGetCorners(user->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL); CHKERRQ(info);
  info = DAGetGhostCorners(user->da,&gxs,&gys,PETSC_NULL,&gxm,&gym,PETSC_NULL); CHKERRQ(info);

  /* Scatter ghost points to local vector */
  info = DAGlobalToLocalBegin(user->da,X,INSERT_VALUES,localX); CHKERRQ(info);
  info = DAGlobalToLocalEnd(user->da,X,INSERT_VALUES,localX); CHKERRQ(info);

  /* Get pointers to vector data */
  info = DAVecGetArray(user->da,localX,(void**)&x);
  info = DAVecGetArray(user->da,G,(void**)&g);

  /* Compute function and gradient over the locally owned part of the mesh */
  for (j=ys; j<ys+ym; j++){
    for (i=xs; i< xs+xm; i++){
      
      xc = x[j][i];
      xlt=xrb=xl=xr=xb=xt=xc;
      
      if (i==0){ /* left side */
        xl= user->left[j-ys+1];
        xlt = user->left[j-ys+2];
      } else {
        xl = x[j][i-1];
      }

      if (j==0){ /* bottom side */
        xb=user->bottom[i-xs+1];
        xrb =user->bottom[i-xs+2];
      } else {
        xb = x[j-1][i];
      }
      
      if (i+1 == gxs+gxm){ /* right side */
        xr=user->right[j-ys+1];
        xrb = user->right[j-ys];
      } else {
        xr = x[j][i+1];
      }

      if (j+1==gys+gym){ /* top side */
        xt=user->top[i-xs+1];
        xlt = user->top[i-xs];
      }else {
        xt = x[j+1][i];
      }

      if (i>gxs && j+1<gys+gym){
        xlt = x[j+1][i-1];
      }
      if (j>gys && i+1<gxs+gxm){
        xrb = x[j-1][i+1];
      }

      d1 = (xc-xl);
      d2 = (xc-xr);
      d3 = (xc-xt);
      d4 = (xc-xb);
      d5 = (xr-xrb);
      d6 = (xrb-xb);
      d7 = (xlt-xl);
      d8 = (xt-xlt);
      
      df1dxc = d1*hydhx;
      df2dxc = ( d1*hydhx + d4*hxdhy );
      df3dxc = d3*hxdhy;
      df4dxc = ( d2*hydhx + d3*hxdhy );
      df5dxc = d2*hydhx;
      df6dxc = d4*hxdhy;

      d1 *= rhx;
      d2 *= rhx;
      d3 *= rhy;
      d4 *= rhy;
      d5 *= rhy;
      d6 *= rhx;
      d7 *= rhy;
      d8 *= rhx;

      f1 = sqrt( 1.0 + d1*d1 + d7*d7);
      f2 = sqrt( 1.0 + d1*d1 + d4*d4);
      f3 = sqrt( 1.0 + d3*d3 + d8*d8);
      f4 = sqrt( 1.0 + d3*d3 + d2*d2);
      f5 = sqrt( 1.0 + d2*d2 + d5*d5);
      f6 = sqrt( 1.0 + d4*d4 + d6*d6);
      
      f2 = sqrt( 1.0 + d1*d1 + d4*d4);
      f4 = sqrt( 1.0 + d3*d3 + d2*d2);

      ft = ft + (f2 + f4);

      df1dxc /= f1;
      df2dxc /= f2;
      df3dxc /= f3;
      df4dxc /= f4;
      df5dxc /= f5;
      df6dxc /= f6;

      g[j][i] = (df1dxc+df2dxc+df3dxc+df4dxc+df5dxc+df6dxc ) * 0.5;
      
    }
  }

  /* Compute triangular areas along the border of the domain. */
  if (xs==0){ /* left side */
    for (j=ys; j<ys+ym; j++){
      d3=(user->left[j-ys+1] - user->left[j-ys+2])*rhy;
      d2=(user->left[j-ys+1] - x[j][0]) *rhx;
      ft = ft+sqrt( 1.0 + d3*d3 + d2*d2);
    }
  }
  if (ys==0){ /* bottom side */
    for (i=xs; i<xs+xm; i++){
      d2=(user->bottom[i+1-xs]-user->bottom[i-xs+2])*rhx;
      d3=(user->bottom[i-xs+1]-x[0][i])*rhy;
      ft = ft+sqrt( 1.0 + d3*d3 + d2*d2);
    }
  }

  if (xs+xm==mx){ /* right side */
    for (j=ys; j< ys+ym; j++){
      d1=(x[j][mx-1] - user->right[j-ys+1])*rhx;
      d4=(user->right[j-ys]-user->right[j-ys+1])*rhy;
      ft = ft+sqrt( 1.0 + d1*d1 + d4*d4);
    }
  }
  if (ys+ym==my){ /* top side */
    for (i=xs; i<xs+xm; i++){
      d1=(x[my-1][i] - user->top[i-xs+1])*rhy;
      d4=(user->top[i-xs+1] - user->top[i-xs])*rhx;
      ft = ft+sqrt( 1.0 + d1*d1 + d4*d4);
    }
  }

  if (ys==0 && xs==0){
    d1=(user->left[0]-user->left[1])/hy;
    d2=(user->bottom[0]-user->bottom[1])*rhx;
    ft +=sqrt( 1.0 + d1*d1 + d2*d2);
  }
  if (ys+ym == my && xs+xm == mx){
    d1=(user->right[ym+1] - user->right[ym])*rhy;
    d2=(user->top[xm+1] - user->top[xm])*rhx;
    ft +=sqrt( 1.0 + d1*d1 + d2*d2);
  }

  ft=ft*area;
  info = MPI_Allreduce(&ft,fcn,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);CHKERRQ(info);

  /* Restore vectors */
  info = DAVecRestoreArray(user->da,localX,(void**)&x);
  info = DAVecRestoreArray(user->da,G,(void**)&g);

  /* Scatter values to global vector */
  info = DARestoreLocalVector(user->da,&localX); CHKERRQ(info);

  info = PetscLogFlops(67*xm*ym); CHKERRQ(info);

  return 0;
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormHessian"
/*
   FormHessian - Evaluates Hessian matrix.

   Input Parameters:
.  taoapp  - the TAO_APPLICATION context
.  x    - input vector
.  ptr  - optional user-defined context, as set by TaoSetHessian()

   Output Parameters:
.  H    - Hessian matrix
.  Hpre - optionally different preconditioning matrix
.  flg  - flag indicating matrix structure

*/
int FormHessian(TAO_APPLICATION taoapp,Vec X,Mat *H, Mat *Hpre, MatStructure *flg, void *ptr)
{ 
  int    info;
  AppCtx *user = (AppCtx *) ptr;

  /* Evaluate the Hessian entries*/
  info = QuadraticH(user,X,*H); CHKERRQ(info);


  /* Indicate that this matrix has the same sparsity pattern during
     successive iterations; setting this flag can save significant work
     in computing the preconditioner for some methods. */
  *flg=SAME_NONZERO_PATTERN;

  return 0;
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "QuadraticH"
/*
   QuadraticH - Evaluates Hessian matrix.

   Input Parameters:
.  user - user-defined context, as set by TaoSetHessian()
.  X    - input vector

   Output Parameter:
.  H    - Hessian matrix
*/
int QuadraticH(AppCtx *user, Vec X, Mat Hessian)
{
  int    info;
  PetscInt i,j,k;
  PetscInt mx=user->mx, my=user->my;
  PetscInt xs,xm,gxs,gxm,ys,ym,gys,gym;
  double hx=1.0/(mx+1), hy=1.0/(my+1), hydhx=hy/hx, hxdhy=hx/hy;
  double f1,f2,f3,f4,f5,f6,d1,d2,d3,d4,d5,d6,d7,d8,xc,xl,xr,xt,xb,xlt,xrb;
  double hl,hr,ht,hb,hc,htl,hbr;
  PetscScalar **x, v[7];
  MatStencil col[7],row;
  Vec    localX;
  PetscTruth assembled;

  /* Get local mesh boundaries */
  info = DAGetLocalVector(user->da,&localX);CHKERRQ(info);

  info = DAGetCorners(user->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL); CHKERRQ(info);
  info = DAGetGhostCorners(user->da,&gxs,&gys,PETSC_NULL,&gxm,&gym,PETSC_NULL); CHKERRQ(info);

  /* Scatter ghost points to local vector */
  info = DAGlobalToLocalBegin(user->da,X,INSERT_VALUES,localX); CHKERRQ(info);
  info = DAGlobalToLocalEnd(user->da,X,INSERT_VALUES,localX); CHKERRQ(info);

  /* Get pointers to vector data */
  info = DAVecGetArray(user->da,localX,(void**)&x);

  /* Initialize matrix entries to zero */
  info = MatAssembled(Hessian,&assembled); CHKERRQ(info);
  if (assembled){info = MatZeroEntries(Hessian);  CHKERRQ(info);}


  /* Set various matrix options */
  info = MatSetOption(Hessian,MAT_IGNORE_OFF_PROC_ENTRIES,PETSC_TRUE); CHKERRQ(info);

  /* Compute Hessian over the locally owned part of the mesh */

  for (j=ys; j<ys+ym; j++){
      
    for (i=xs; i< xs+xm; i++){

      xc = x[j][i];
      xlt=xrb=xl=xr=xb=xt=xc;

      /* Left side */
      if (i==0){
        xl  = user->left[j-ys+1];
        xlt = user->left[j-ys+2];
      } else {
        xl  = x[j][i-1];
      }
      
      if (j==0){
        xb  = user->bottom[i-xs+1];
        xrb = user->bottom[i-xs+2];
      } else {
        xb  = x[j-1][i];
      }
      
      if (i+1 == mx){
        xr  = user->right[j-ys+1];
        xrb = user->right[j-ys];
      } else {
        xr  = x[j][i+1];
      }

      if (j+1==my){
        xt  = user->top[i-xs+1];
        xlt = user->top[i-xs];
      }else {
        xt  = x[j+1][i];
      }

      if (i>0 && j+1<my){
        xlt = x[j+1][i-1];
      }
      if (j>0 && i+1<mx){
        xrb = x[j-1][i+1];
      }


      d1 = (xc-xl)/hx;
      d2 = (xc-xr)/hx;
      d3 = (xc-xt)/hy;
      d4 = (xc-xb)/hy;
      d5 = (xrb-xr)/hy;
      d6 = (xrb-xb)/hx;
      d7 = (xlt-xl)/hy;
      d8 = (xlt-xt)/hx;
      
      f1 = sqrt( 1.0 + d1*d1 + d7*d7);
      f2 = sqrt( 1.0 + d1*d1 + d4*d4);
      f3 = sqrt( 1.0 + d3*d3 + d8*d8);
      f4 = sqrt( 1.0 + d3*d3 + d2*d2);
      f5 = sqrt( 1.0 + d2*d2 + d5*d5);
      f6 = sqrt( 1.0 + d4*d4 + d6*d6);


      hl = (-hydhx*(1.0+d7*d7)+d1*d7)/(f1*f1*f1)+
	(-hydhx*(1.0+d4*d4)+d1*d4)/(f2*f2*f2);
      hr = (-hydhx*(1.0+d5*d5)+d2*d5)/(f5*f5*f5)+
	(-hydhx*(1.0+d3*d3)+d2*d3)/(f4*f4*f4);
      ht = (-hxdhy*(1.0+d8*d8)+d3*d8)/(f3*f3*f3)+
	(-hxdhy*(1.0+d2*d2)+d2*d3)/(f4*f4*f4);
      hb = (-hxdhy*(1.0+d6*d6)+d4*d6)/(f6*f6*f6)+
	(-hxdhy*(1.0+d1*d1)+d1*d4)/(f2*f2*f2);

      hbr = -d2*d5/(f5*f5*f5) - d4*d6/(f6*f6*f6);
      htl = -d1*d7/(f1*f1*f1) - d3*d8/(f3*f3*f3);

      hc = hydhx*(1.0+d7*d7)/(f1*f1*f1) + hxdhy*(1.0+d8*d8)/(f3*f3*f3) +
	hydhx*(1.0+d5*d5)/(f5*f5*f5) + hxdhy*(1.0+d6*d6)/(f6*f6*f6) +
	(hxdhy*(1.0+d1*d1)+hydhx*(1.0+d4*d4)-2*d1*d4)/(f2*f2*f2) +
	(hxdhy*(1.0+d2*d2)+hydhx*(1.0+d3*d3)-2*d2*d3)/(f4*f4*f4);

      hl/=2.0; hr/=2.0; ht/=2.0; hb/=2.0; hbr/=2.0; htl/=2.0;  hc/=2.0; 

      row.j = j; row.i = i;
      k=0;
      if (j>0){ 
	v[k]=hb;
	col[k].j = j - 1; col[k].i = i;
	k++;
      }
      
      if (j>0 && i < mx -1){
	v[k]=hbr;
	col[k].j = j - 1; col[k].i = i+1;
	k++;
      }
      
      if (i>0){
	v[k]= hl;
	col[k].j = j; col[k].i = i-1;
	k++;
      }
      
      v[k]= hc;
      col[k].j = j; col[k].i = i;
      k++;
      
      if (i < mx-1 ){
	v[k]= hr; 
	col[k].j = j; col[k].i = i+1;
	k++;
      }
      
      if (i>0 && j < my-1 ){
	v[k]= htl;
	col[k].j = j+1; col[k].i = i-1;
	k++;
      }
      
      if (j < my-1 ){
	v[k]= ht; 
	col[k].j = j+1; col[k].i = i;
	k++;
      }
      
      /* 
	 Set matrix values using local numbering, which was defined
	 earlier, in the main routine.
      */
      info = MatSetValuesStencil(Hessian,1,&row,k,col,v,INSERT_VALUES);
      CHKERRQ(info);
      
    }
  }
  
  /* Restore vectors */
  info = DAVecRestoreArray(user->da,localX,(void**)&x);

  info = DARestoreLocalVector(user->da,&localX); CHKERRQ(info);

  /* Assemble the matrix */
  info = MatAssemblyBegin(Hessian,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  info = MatAssemblyEnd(Hessian,MAT_FINAL_ASSEMBLY); CHKERRQ(info);

  info = PetscLogFlops(199*xm*ym); CHKERRQ(info);
  return 0;
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "MSA_BoundaryConditions"
/* 
   MSA_BoundaryConditions -  Calculates the boundary conditions for
   the region.

   Input Parameter:
.  user - user-defined application context

   Output Parameter:
.  user - user-defined application context
*/
static int MSA_BoundaryConditions(AppCtx * user)
{
  int        info;
  PetscInt     i,j,k,limit=0,maxits=5;
  PetscInt     xs,ys,xm,ym,gxs,gys,gxm,gym;
  PetscInt     mx=user->mx,my=user->my;
  PetscInt     bsize=0, lsize=0, tsize=0, rsize=0;
  double     one=1.0, two=2.0, three=3.0, tol=1e-10;
  double     fnorm,det,hx,hy,xt=0,yt=0;
  double     u1,u2,nf1,nf2,njac11,njac12,njac21,njac22;
  double     b=-0.5, t=0.5, l=-0.5, r=0.5;
  double     *boundary;
  PetscTruth   flg;

  /* Get local mesh boundaries */
  info = DAGetCorners(user->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL); CHKERRQ(info);
  info = DAGetGhostCorners(user->da,&gxs,&gys,PETSC_NULL,&gxm,&gym,PETSC_NULL); CHKERRQ(info);

  bsize=xm+2;
  lsize=ym+2;
  rsize=ym+2;
  tsize=xm+2;

  info = PetscMalloc(bsize*sizeof(double),&user->bottom); CHKERRQ(info);
  info = PetscMalloc(tsize*sizeof(double),&user->top); CHKERRQ(info);
  info = PetscMalloc(lsize*sizeof(double),&user->left); CHKERRQ(info);
  info = PetscMalloc(rsize*sizeof(double),&user->right); CHKERRQ(info);

  hx= (r-l)/(mx+1); hy=(t-b)/(my+1);

  for (j=0; j<4; j++){
    if (j==0){
      yt=b;
      xt=l+hx*xs;
      limit=bsize;
      boundary=user->bottom;
    } else if (j==1){
      yt=t;
      xt=l+hx*xs;
      limit=tsize;
      boundary=user->top;
    } else if (j==2){
      yt=b+hy*ys;
      xt=l;
      limit=lsize;
      boundary=user->left;
    } else { //if (j==3)
      yt=b+hy*ys;
      xt=r;
      limit=rsize;
      boundary=user->right;
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

      boundary[i]=u1*u1-u2*u2;
      if (j==0 || j==1) {
	xt=xt+hx;
      } else { // if (j==2 || j==3)
	yt=yt+hy;
      }
      
    }

  }

  /* Scale the boundary if desired */
  if (1==1){ 
    PetscReal scl = 1.0;

    info = PetscOptionsGetReal(PETSC_NULL,"-bottom",&scl,&flg); 
    CHKERRQ(info);
    if (flg){
      for (i=0;i<bsize;i++) user->bottom[i]*=scl;
    }

    info = PetscOptionsGetReal(PETSC_NULL,"-top",&scl,&flg); 
    CHKERRQ(info);
    if (flg){
      for (i=0;i<tsize;i++) user->top[i]*=scl;
    }

    info = PetscOptionsGetReal(PETSC_NULL,"-right",&scl,&flg); 
    CHKERRQ(info);
    if (flg){
      for (i=0;i<rsize;i++) user->right[i]*=scl;
    }

    info = PetscOptionsGetReal(PETSC_NULL,"-left",&scl,&flg); 
    CHKERRQ(info);
    if (flg){
      for (i=0;i<lsize;i++) user->left[i]*=scl;
    }
  }
  
  return 0;
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "MSA_InitialPoint"
/*
   MSA_InitialPoint - Calculates the initial guess in one of three ways. 

   Input Parameters:
.  user - user-defined application context
.  X - vector for initial guess

   Output Parameters:
.  X - newly computed initial guess
*/
static int MSA_InitialPoint(AppCtx * user, Vec X)
{
  int      info;
  PetscInt   start2=-1,i,j;
  PetscReal   start1=0;
  PetscTruth flg1,flg2;

  info = PetscOptionsGetReal(PETSC_NULL,"-start",&start1,&flg1); CHKERRQ(info);
  info = PetscOptionsGetInt(PETSC_NULL,"-random",&start2,&flg2); CHKERRQ(info);

  if (flg1){ /* The zero vector is reasonable */
 
    info = VecSet(X, start1); CHKERRQ(info);

  } else if (flg2 && start2>0){ /* Try a random start between -0.5 and 0.5 */

    PetscRandom rctx;  PetscScalar np5=-0.5;

    info = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);
    CHKERRQ(info);
    for (i=0; i<start2; i++){
      info = VecSetRandom(X, rctx); CHKERRQ(info);
    }
    info = PetscRandomDestroy(rctx); CHKERRQ(info);
    info = VecShift(X, np5); CHKERRQ(info);

  } else { /* Take an average of the boundary conditions */

    PetscInt xs,xm,ys,ym;
    PetscInt mx=user->mx,my=user->my;
    PetscScalar **x;
    
    /* Get local mesh boundaries */
    info = DAGetCorners(user->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL); CHKERRQ(info);
    
    /* Get pointers to vector data */
    info = DAVecGetArray(user->da,X,(void**)&x);

    /* Perform local computations */    
    for (j=ys; j<ys+ym; j++){
      for (i=xs; i< xs+xm; i++){
	x[j][i] = ( ((j+1)*user->bottom[i-xs+1]+(my-j+1)*user->top[i-xs+1])/(my+2)+
		   ((i+1)*user->left[j-ys+1]+(mx-i+1)*user->right[j-ys+1])/(mx+2))/2.0; 
      }
    }
    
    /* Restore vectors */
    info = DAVecRestoreArray(user->da,X,(void**)&x);  CHKERRQ(info);

    info = PetscLogFlops(9*xm*ym); CHKERRQ(info);
    
  }
  return 0;
}

/*-----------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "My_Monitor"
int My_Monitor(TAO_APPLICATION minsurfapp, void *ctx){
  int info;
  Vec X;

  info = TaoAppGetSolutionVec(minsurfapp,&X); CHKERRQ(info);
  info = VecView(X,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(info);
  return 0;
}
